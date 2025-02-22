/*
 *
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 * Written by Slaven Peles <peles2@llnl.gov>.
 * LLNL-CODE-718378.
 * All rights reserved.
 *
 * This file is part of GridKit™. For details, see github.com/LLNL/GridKit
 * Please also read the LICENSE file.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * - Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the disclaimer below.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the disclaimer (as noted below) in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the LLNS/LLNL nor the names of its contributors may
 *   be used to endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
 * SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Lawrence Livermore National Laboratory is operated by Lawrence Livermore
 * National Security, LLC, for the U.S. Department of Energy, National
 * Nuclear Security Administration under Contract DE-AC52-07NA27344.
 *
 * This document was prepared as an account of work sponsored by an agency
 * of the United States government. Neither the United States government nor
 * Lawrence Livermore National Security, LLC, nor any of their employees
 * makes any warranty, expressed or implied, or assumes any legal liability
 * or responsibility for the accuracy, completeness, or usefulness of any
 * information, apparatus, product, or process disclosed, or represents that
 * its use would not infringe privately owned rights. Reference herein to
 * any specific commercial product, process, or service by trade name,
 * trademark, manufacturer, or otherwise does not necessarily constitute or
 * imply its endorsement, recommendation, or favoring by the United States
 * government or Lawrence Livermore National Security, LLC. The views and
 * opinions of authors expressed herein do not necessarily state or reflect
 * those of the United States government or Lawrence Livermore National
 * Security, LLC, and shall not be used for advertising or product
 * endorsement purposes.
 *
 */


#include <iostream>
#include <iomanip>

#include <idas/idas.h>
#include <idas/idas_ls.h>

#include "Evaluator.hpp"
#include "Ida.hpp"


namespace AnalysisManager
{

namespace Sundials
{

    template <class ScalarT, typename IdxT>
    Ida<ScalarT, IdxT>::Ida(GridKit::Model::Evaluator<ScalarT, IdxT>* model) : DynamicSolver<ScalarT, IdxT>(model)
    {
        int retval = 0;

        // Create the SUNDIALS context that all SUNDIALS objects require
        retval = SUNContext_Create(SUN_COMM_NULL, &context_);
        checkOutput(retval, "SUNContext");
        solver_ = IDACreate(context_);
        tag_ = NULL;
    }

    /**
     * @brief Destroy the Ida< Scalar T,  Idx T>:: Ida object
     * 
     * @note if sysmodel is freed before this will fail. May want something agnostic to this
     * 
     * @tparam ScalarT 
     * @tparam IdxT 
     */
    template <class ScalarT, typename IdxT>
    Ida<ScalarT, IdxT>::~Ida()
    {
        N_VDestroy(yy_);
        N_VDestroy(yp_);
        N_VDestroy(yy0_);
        N_VDestroy(yp0_);
        if (model_->hasJacobian())
        {
            SUNLinSolFree_KLU(linearSolver_);
            SUNMatDestroy_Sparse(JacobianMat_);
        }
        else
        {
            SUNLinSolFree_Dense(linearSolver_);
            SUNMatDestroy_Dense(JacobianMat_);
        }
        ///@todo this free is needed but on geninfbus this seg faults
        //IDAFree(&solver_);
        SUNContext_Free(&context_);
        
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureSimulation()
    {
        int retval = 0;

        // Allocate solution vectors
        yy_ = N_VNew_Serial(model_->size(), context_);
        checkAllocation((void*) yy_, "N_VNew_Serial");
        yp_ = N_VClone(yy_);
        checkAllocation((void*) yp_, "N_VClone");

        //get intial conditions
        this->getDefaultInitialCondition();

        // Create vectors to store restart initial condition
        yy0_ = N_VClone(yy_);
        checkAllocation((void*) yy0_, "N_VClone");
        yp0_ = N_VClone(yp_);
        checkAllocation((void*) yp0_, "N_VClone");

        // Dummy initial time; will be overridden.
        const sunrealtype t0 = SUN_RCONST(0.0);

        // Allocate and initialize IDA workspace
        retval = IDAInit(solver_, this->Residual, t0, yy_, yp_);
        checkOutput(retval, "IDAInit");

        // Set pointer to model data
        retval = IDASetUserData(solver_, model_);
        checkOutput(retval, "IDASetUserData");

        // Set tolerances
        sunrealtype rel_tol;
        sunrealtype abs_tol;

        model_->setTolerances(rel_tol, abs_tol); ///< \todo Function name should be "getTolerances"!
        retval = IDASStolerances(solver_, rel_tol, abs_tol);
        checkOutput(retval, "IDASStolerances");
        
        IdxT msa;
        model_->setMaxSteps(msa);

        /// \todo Need to set max number of steps based on user input!
        retval = IDASetMaxNumSteps(solver_, msa);
        checkOutput(retval, "IDASetMaxNumSteps");

        // Tag differential variables
        std::vector<bool>& tag = model_->tag();
        if (static_cast<IdxT>(tag.size()) == model_->size())
        {
            tag_ = N_VClone(yy_);
            checkAllocation((void*) tag_, "N_VClone");
            model_->tagDifferentiable();
            copyVec(tag, tag_);

            retval = IDASetId(solver_, tag_);
            checkOutput(retval, "IDASetId");
            retval = IDASetSuppressAlg(solver_, SUNTRUE);
            checkOutput(retval, "IDASetSuppressAlg");
        }

        // Set up linear solver
        this->configureLinearSolver();

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolver()
    {
        int retval = 0;
        if (model_->hasJacobian())
        {
            JacobianMat_ = SUNSparseMatrix(model_->size(), model_->size(), model_->size() * model_->size(), CSR_MAT, context_);
            checkAllocation((void*) JacobianMat_, "SUNSparseMatrix");

            linearSolver_ = SUNLinSol_KLU(yy_, JacobianMat_, context_);
            checkAllocation((void*) linearSolver_, "SUNLinSol_KLU");

            retval = IDASetLinearSolver(solver_, linearSolver_, JacobianMat_);
            checkOutput(retval, "IDASetLinearSolver");

            retval = IDASetJacFn(solver_, this->Jac);
            checkOutput(retval, "IDASetJacFn");
        }
        else
        {
            JacobianMat_ = SUNDenseMatrix(model_->size(), model_->size(), context_);
            checkAllocation((void*) JacobianMat_, "SUNDenseMatrix");

            linearSolver_ = SUNLinSol_Dense(yy_, JacobianMat_, context_);
            checkAllocation((void*) linearSolver_, "SUNLinSol_Dense");

            retval = IDASetLinearSolver(solver_, linearSolver_, JacobianMat_);
            checkOutput(retval, "IDASetLinearSolver");

        }

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::getDefaultInitialCondition()
    {
        model_->initialize();

        copyVec(model_->y(), yy_);
        copyVec(model_->yp(), yp_);

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::setIntegrationTime(real_type t_init, real_type t_final, int nout)
    {
        t_init_  = t_init;
        t_final_ = t_final;
        nout_    = nout;
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeSimulation(real_type t0, bool findConsistent)
    {
        int retval = 0;

        // Need to reinitialize IDA to set to get correct initial conditions
        retval = IDAReInit(solver_, t0, yy_, yp_);
        checkOutput(retval, "IDAReInit");

        // Find a consistent set of initial conditions for DAE
        if (findConsistent)
        {
            int initType = IDA_Y_INIT;

            if (tag_)
                initType = IDA_YA_YDP_INIT;

            retval = IDACalcIC(solver_, initType, 0.1);
            checkOutput(retval, "IDACalcIC");
            
            copyVec(yy_, model_->y());
            copyVec(yp_, model_->yp());
        }

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runSimulationFixed(real_type t0, real_type dt, 
        real_type tmax, FILE *f)
    {
        int retval = 0;
        int iout = 0;
        real_type t, tret;
        
        //printOutputF(t0, 0, f);
        for (t = t0+dt; t < tmax; t += dt)
        {
            retval = IDASolve(solver_, t, &tret, yy_, yp_, IDA_NORMAL);
            checkOutput(retval, "IDASolve");
            printOutputF(t, retval, f); 

            if (retval != IDA_SUCCESS)
            {
                printf("IDA Failure! %dn", retval);
                break;
            }
        }

        model_->updateTime(t, 0.0);
        copyVec(yy_, model_->y());
        copyVec(yp_, model_->yp());

        return retval;
    }


    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runSimulation(real_type tf, int nout)
    {
        int retval = 0;
        int iout = 0;
        real_type tret;
        real_type dt = tf / static_cast<real_type>(nout);
        real_type tout = dt;

        /* In loop, call IDASolve, print results, and test for error.
         *     Break out of loop when NOUT preset output times have been reached. */
        printOutput(0.0);
        while(nout > iout)
        {
            retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
            checkOutput(retval, "IDASolve");
            printOutput(tout); 

            if (retval == IDA_SUCCESS)
            {
                ++iout;
                tout += dt;
            }
        }

        //Final copy out. No gaurentee last residual evaluation is final step.
        model_->updateTime(tf, 0.0);
        copyVec(yy_, model_->y());
        copyVec(yp_, model_->yp());

        //std::cout << "\n";
        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteSimulation()
    {
        IDAFree(&solver_);
        SUNLinSolFree(linearSolver_);
        N_VDestroy(yy_);
        N_VDestroy(yp_);
        return 0;
    }


    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureQuadrature()
    {
        int retval = 0;

        // Create and initialize quadratures
        q_ = N_VNew_Serial(model_->sizeQuadrature(), context_);
        checkAllocation((void*) q_, "N_VNew_Serial");

        // Set integrand function and allocate quadrature workspace
        retval = IDAQuadInit(solver_, this->Integrand, q_);
        checkOutput(retval, "IDAQuadInit");

        // Set tolerances and error control for quadratures
        real_type rel_tol, abs_tol;
        model_->setTolerances(rel_tol, abs_tol);

        // Set tolerances for quadrature stricter than for integration
        retval = IDAQuadSStolerances(solver_, rel_tol*0.1, abs_tol*0.1);
        checkOutput(retval, "IDAQuadSStolerances");

        // Include quadrature in eror checking
        retval = IDASetQuadErrCon(solver_, SUNTRUE);
        checkOutput(retval, "IDASetQuadErrCon");

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeQuadrature()
    {
        int retval = 0;

        // Set all quadratures to zero
        N_VConst(SUN_RCONST(0.0), q_);

        // Initialize quadratures
        retval = IDAQuadReInit(solver_, q_);
        checkOutput(retval, "IDAQuadInit");

        return retval;
    }


    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runSimulationQuadrature(real_type tf, int nout)
    {
        int retval = 0;
        real_type tret;

        //std::cout << "Forward integration for initial value problem ... \n";

        real_type dt = tf / static_cast<real_type>(nout);
        real_type tout = dt;
        //printOutput(0.0);
        //printSpecial(0.0, yy_);
        for(int i = 0; i < nout; ++i)
        {
            retval = IDASolve(solver_, tout, &tret, yy_, yp_, IDA_NORMAL);
            checkOutput(retval, "IDASolve");
            //printSpecial(tout, yy_);
            //printOutput(tout); 

            if (retval == IDA_SUCCESS)
            {
                tout += dt;
            }

            retval = IDAGetQuad(solver_, &tret, q_);
            checkOutput(retval, "IDAGetQuad");
        }
        
        //Final copy out. No gaurentee last residual evaluation is final step.
        model_->updateTime(tf, 0.0);
        copyVec(yy_, model_->y());
        copyVec(yp_, model_->yp());

        return retval;
    }


    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteQuadrature()
    {
        IDAQuadFree(solver_);
        N_VDestroy(q_);

        return 0;
    }


    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureAdjoint()
    {
        // Allocate adjoint vector, derivatives and quadrature
      yyB_ = N_VNew_Serial(model_->size(), context_);
        checkAllocation((void*) yyB_, "N_VNew_Serial");

        ypB_ = N_VClone(yyB_);
        checkAllocation((void*) ypB_, "N_VClone");

        qB_ = N_VNew_Serial(model_->sizeParams(), context_);
        checkAllocation((void*) qB_, "N_VNew_Serial");

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeAdjoint(IdxT steps)
    {
        int retval = 0;

        // Create adjoint workspace
        retval = IDAAdjInit(solver_, steps, IDA_HERMITE);
        checkOutput(retval, "IDAAdjInit");

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::initializeBackwardSimulation(real_type tf)
    {
        int retval = 0;
        sunrealtype rel_tol;
        sunrealtype abs_tol;

        model_->initializeAdjoint();

        copyVec(model_->yB(),  yyB_);
        copyVec(model_->ypB(), ypB_);
        N_VConst(0.0, qB_);

        retval = IDACreateB(solver_, &backwardID_);
        checkOutput(retval, "IDACreateB");

        // IDAInitB must be called after forward simulation run.
        retval = IDAInitB(solver_, backwardID_, this->adjointResidual, tf, yyB_, ypB_);
        checkOutput(retval, "IDAInitB");

        model_->setTolerances(rel_tol, abs_tol);
        retval = IDASStolerancesB(solver_, backwardID_, rel_tol, abs_tol);
        checkOutput(retval, "IDASStolerancesB");

        retval = IDASetUserDataB(solver_, backwardID_, model_);
        checkOutput(retval, "IDASetUserDataB");

        /// \todo Need to set max number of steps based on user input!
        retval = IDASetMaxNumStepsB(solver_, backwardID_, 2000);
        checkOutput(retval, "IDASetMaxNumSteps");

        // Set up linear solver
        JacobianMatB_ = SUNDenseMatrix(model_->size(), model_->size(), context_);
        checkAllocation((void*) JacobianMatB_, "SUNDenseMatrix");

        linearSolverB_ = SUNLinSol_Dense(yyB_, JacobianMatB_, context_);
        checkAllocation((void*) linearSolverB_, "SUNLinSol_Dense");

        retval = IDASetLinearSolverB(solver_, backwardID_, linearSolverB_, JacobianMatB_);
        checkOutput(retval, "IDASetLinearSolverB");


        // Also reinitialize quadratures.
        retval = IDAQuadInitB(solver_, backwardID_, this->adjointIntegrand, qB_);
        checkOutput(retval, "IDAQuadInitB");

        //retval = IDAQuadSStolerancesB(solver_, backwardID_, rel_tol*1.1, abs_tol*1.1);
        retval = IDAQuadSStolerancesB(solver_, backwardID_, rel_tol*0.1, abs_tol*0.1);
        checkOutput(retval, "IDAQuadSStolerancesB");

        // Include quadratures in error control
        retval = IDASetQuadErrConB(solver_, backwardID_, SUNTRUE);
        checkOutput(retval, "IDASetQuadErrConB");


        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::configureLinearSolverBackward()
    {
        int retval = 0;

        // Create Jacobian matrix
        JacobianMatB_ = SUNDenseMatrix(model_->size(), model_->size(), context_);
        checkAllocation((void*) JacobianMatB_, "SUNDenseMatrix");

        // Create linear solver
        linearSolverB_ = SUNLinSol_Dense(yyB_, JacobianMatB_, context_);
        checkAllocation((void*) linearSolverB_, "SUNLinSol_Dense");

        // Attach linear solver to IDA
        retval = IDASetLinearSolverB(solver_, backwardID_, linearSolverB_, JacobianMatB_);
        checkOutput(retval, "IDASetLinearSolverB");

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runForwardSimulation(real_type tf, int nout)
    {
        int retval = 0;
        int ncheck;
        real_type time;

        //std::cout << "Forward integration for adjoint analysis ... \n";

        real_type dt = tf / static_cast<real_type>(nout);
        real_type tout = dt;
        for(int i = 0; i < nout; ++i)
        {
            retval = IDASolveF(solver_, tout, &time, yy_, yp_, IDA_NORMAL, &ncheck);
            checkOutput(retval, "IDASolveF");

            if (retval == IDA_SUCCESS)
            {
                tout += dt;
            }

            retval = IDAGetQuad(solver_, &time, q_);
            checkOutput(retval, "IDASolve");
        }
        
        //Final copy out. No gaurentee last residual evaluation is final step.
        model_->updateTime(tf, 0.0);
        copyVec(yy_, model_->y());
        copyVec(yp_, model_->yp());

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::runBackwardSimulation(real_type t_init)
    {
        int retval = 0;
        long int nstB;
        real_type time;

        //std::cout << "Backward integration for adjoint analysis ... ";

        retval = IDASolveB(solver_, t_init, IDA_NORMAL);
        checkOutput(retval, "IDASolveB");

        IDAGetNumSteps(IDAGetAdjIDABmem(solver_, backwardID_), &nstB);
        //std::cout << "done ( nst = " << nstB << " )\n";

        retval = IDAGetB(solver_, backwardID_, &time, yyB_, ypB_);
        checkOutput(retval, "IDAGetB");

        //Copy back into model
        copyVec(yyB_, model_->yB());
        copyVec(ypB_, model_->ypB());

        retval = IDAGetQuadB(solver_, backwardID_, &time, qB_);
        checkOutput(retval, "IDAGetQuadB");

        return retval;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::deleteAdjoint()
    {
        IDAAdjFree(solver_);
        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Residual(sunrealtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
    {
        GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

        model->updateTime(tres, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());

        model->evaluateResidual();
        const std::vector<ScalarT>& f = model->getResidual();
        copyVec(f, rr);

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Jac(sunrealtype t, sunrealtype cj, N_Vector yy, N_Vector yp, N_Vector resvec, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
    {

        GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);


        model->updateTime(t, cj);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());

        model->evaluateJacobian();
        COO_Matrix<ScalarT, IdxT> Jac = model->getJacobian();
        
        //Get reference to the jacobian entries
        std::tuple<std::vector<IdxT>&, std::vector<IdxT>&, std::vector<ScalarT>&> tpm = Jac.getEntries();
        const auto [r, c, val] = tpm;

        //get the CSR row pointers from COO matrix
        std::vector<IdxT> csrrowdata = Jac.getCSRRowData();

        SUNMatZero(J);

        //Set row pointers
        sunindextype *rowptrs = SUNSparseMatrix_IndexPointers(J);
        for (unsigned int i = 0; i < csrrowdata.size() ; i++)
        {
            rowptrs[i] = csrrowdata[i];
        }

        sunindextype *colvals = SUNSparseMatrix_IndexValues(J);
        sunrealtype  *data = SUNSparseMatrix_Data(J);
        //Copy data from model jac to sundials
        for (unsigned int i = 0; i < c.size(); i++ )
        {
            colvals[i] = c[i];
            data[i] = val[i];
        }

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::Integrand(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector rhsQ, void *user_data)
    {
        GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

        model->updateTime(tt, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());

        model->evaluateIntegrand();
        const std::vector<ScalarT>& g = model->getIntegrand();
        copyVec(g, rhsQ);

        return 0;
    }

    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::adjointResidual(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rrB, void *user_data)
    {
        GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

        model->updateTime(tt, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());
        copyVec(yyB, model->yB());
        copyVec(ypB, model->ypB());

        model->evaluateAdjointResidual();
        const std::vector<ScalarT>& fB = model->getAdjointResidual();
        copyVec(fB, rrB);

        return 0;
    }


    template <class ScalarT, typename IdxT>
    int Ida<ScalarT, IdxT>::adjointIntegrand(sunrealtype tt, N_Vector yy, N_Vector yp, N_Vector yyB, N_Vector ypB, N_Vector rhsQB, void *user_data)
    {
        GridKit::Model::Evaluator<ScalarT, IdxT>* model = static_cast<GridKit::Model::Evaluator<ScalarT, IdxT>*>(user_data);

        model->updateTime(tt, 0.0);
        copyVec(yy, model->y());
        copyVec(yp, model->yp());
        copyVec(yyB, model->yB());
        copyVec(ypB, model->ypB());

        model->evaluateAdjointIntegrand();
        const std::vector<ScalarT>& gB = model->getAdjointIntegrand();
        copyVec(gB, rhsQB);

        return 0;
    }


    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const N_Vector x, std::vector< ScalarT >& y)
    {
        const ScalarT* xdata = NV_DATA_S(x);
        for(unsigned int i = 0; i < y.size(); ++i)
        {
            y[i] = xdata[i];
        }
    }


    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const std::vector< ScalarT >& x, N_Vector y)
    {
        ScalarT* ydata = NV_DATA_S(y);
        for(unsigned int i = 0; i < x.size(); ++i)
        {
            ydata[i] = x[i];
        }
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::copyVec(const std::vector< bool >& x, N_Vector y)
    {
        ScalarT* ydata = NV_DATA_S(y);
        for(unsigned int i = 0; i < x.size(); ++i)
        {
            if (x[i])
                ydata[i] = 1.0;
            else
                ydata[i] = 0.0;
        }
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printOutputF(sunrealtype t, int res, FILE *f)
    {
        sunrealtype *yval  = N_VGetArrayPointer_Serial(yy_);
        sunrealtype *ypval = N_VGetArrayPointer_Serial(yp_);

        fprintf(f, "%g,%d", t, res);
        for (IdxT i = 0; i < model_->size(); ++i)
        {
            fprintf(f, ",%g", yval[i]);
        }
        for (IdxT i = 0; i < model_->size(); ++i)
        {
            fprintf(f, ",%g", ypval[i]);
        }
        fprintf(f, "\n");
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printOutput(sunrealtype t)
    {
        sunrealtype *yval  = N_VGetArrayPointer_Serial(yy_);
        sunrealtype *ypval = N_VGetArrayPointer_Serial(yp_);

        std::cout << std::setprecision(5) << std::setw(7) << t;
        for (IdxT i = 0; i < model_->size(); ++i)
        {
            std::cout << "," << yval[i];
        }
        for (IdxT i = 0; i < model_->size(); ++i)
        {
            std::cout << "," << ypval[i];
        }
        std::cout << "\n";
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printSpecial(sunrealtype t, N_Vector y)
    {
        sunrealtype *yval = N_VGetArrayPointer_Serial(y);
        IdxT N = static_cast<IdxT>(N_VGetLength_Serial(y));
        std::cout << "{";
        std::cout << std::setprecision(5) << std::setw(7) << t;
        for (IdxT i = 0; i < N; ++i)
        {
            std::cout << ", " << yval[i];
        }
        std::cout << "},\n";
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::printFinalStats()
    {
        int retval = 0;
        void* mem = solver_;
        long int nst;
        long int nre;
        long int nje;
        long int nni;
        long int netf;
        long int ncfn;

        retval = IDAGetNumSteps(mem, &nst);
        checkOutput(retval, "IDAGetNumSteps");
        retval = IDAGetNumResEvals(mem, &nre);
        checkOutput(retval, "IDAGetNumResEvals");
        retval = IDAGetNumJacEvals(mem, &nje);
        checkOutput(retval, "IDAGetNumJacEvals");
        retval = IDAGetNumNonlinSolvIters(mem, &nni);
        checkOutput(retval, "IDAGetNumNonlinSolvIters");
        retval = IDAGetNumErrTestFails(mem, &netf);
        checkOutput(retval, "IDAGetNumErrTestFails");
        retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
        checkOutput(retval, "IDAGetNumNonlinSolvConvFails");

        // std::cout << "\nFinal Run Statistics: \n\n";
        std::cout << "Number of steps                    = " <<  nst  << "\n";
        std::cout << "Number of residual evaluations     = " <<  nre  << "\n";
        //std::cout << "Number of Jacobian evaluations     = " <<  nje  << "\n";
        std::cout << "Number of nonlinear iterations     = " <<  nni  << "\n";
        std::cout << "Number of error test failures      = " <<  netf << "\n";
        std::cout << "Number of nonlinear conv. failures = " <<  ncfn << "\n";
    }


    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::checkAllocation(void* v, const char* functionName)
    {
        if (v == NULL)
        {
            std::cerr << "\nERROR: Function " << functionName << " failed -- returned NULL pointer!\n\n";
            throw SundialsException();
        }
    }

    template <class ScalarT, typename IdxT>
    void Ida<ScalarT, IdxT>::checkOutput(int retval, const char* functionName)
    {
        if (retval < 0)
        {
            std::cerr << "\nERROR: Function " << functionName << " failed with flag " << retval << "!\n\n";
            throw SundialsException();
        }
    }

    // Compiler will prevent building modules with data type incompatible with sunrealtype
    template class Ida<sunrealtype, long int>;
    template class Ida<sunrealtype, int>;
    template class Ida<sunrealtype, size_t>;

} // namespace Sundials
} // namespace AnalysisManager
