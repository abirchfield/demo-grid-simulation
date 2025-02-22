/* Bus Fault Component - Adam Birchfield */
#pragma once

#include <Component.hpp>
#include <BusBase.hpp>

namespace GridKit
{
namespace PhasorDynamics
{
    using ComponentT = Component<double, int>;
    using BaseBusT = BusBase<double, int>;

    class GENROU : public ComponentT
    {
        using ComponentT::size_;
        using ComponentT::nnz_;
        using ComponentT::time_;
        using ComponentT::alpha_;
        using ComponentT::y_;
        using ComponentT::yp_;
        using ComponentT::tag_;
        using ComponentT::f_;
        using ComponentT::g_;
        using ComponentT::yB_;
        using ComponentT::ypB_;
        using ComponentT::fB_;
        using ComponentT::gB_;
        using ComponentT::param_;


    public:

        GENROU(BaseBusT* bus, int unit_id) : bus_(bus), unit_id_(unit_id),
            p0_(0), q0_(0), busID_(0), H_(3), D_(0), Ra_(0), Tdop_(7), 
            Tdopp_(.04), Tqopp_(.05), Tqop_(.75), Xd_(2.1), Xdp_(0.2), 
            Xdpp_(0.18), Xq_(.5), Xqp_(.5), Xqpp_(.18), Xl_(.15),
            S10_(0), S12_(0)
        {
            size_ = 21;
            set_derived_params();
        }

        GENROU(BaseBusT* bus, int unit_id, double p0, double q0, double H, 
            double D, double Ra, double Tdop, double Tdopp, double Tqopp, 
            double Tqop, double Xd, double Xdp, double Xdpp, double Xq, 
            double Xqp, double Xqpp, double Xl, double S10, double S12) : 
            bus_(bus), unit_id_(unit_id), p0_(p0), q0_(q0), busID_(0), H_(H), 
            D_(D), Ra_(Ra), Tdop_(Tdop), Tdopp_(Tdopp), Tqopp_(Tqopp), 
            Tqop_(Tqop), Xd_(Xd), Xdp_(Xdp), Xdpp_(Xdpp), Xq_(Xq), Xqp_(Xqp), 
            Xqpp_(Xqpp), Xl_(Xl), S10_(S10), S12_(S12)
        {
            size_ = 21;
            set_derived_params();
        }

        void set_derived_params() 
        {
            SA_ = 0;
            SB_ = 0;
            if (S12_ != 0)
            {
                double s112 = sqrt(S10_ / S12_);
                SA_ = (1.2*s112 + 1) / (s112 + 1);
                SB_ = (1.2*s112 - 1) / (s112 - 1);
                if (SB_ < SA_) SA_ = SB_;
                SB_ = S12_ / pow(SA_ - 1.2, 2);
            }
            Xd1_ = Xd_ - Xdp_;
            Xd2_ = Xdp_ - Xl_;
            Xd3_ = (Xdp_ - Xdpp_) / (Xd2_ * Xd2_);
            Xd4_ = (Xdp_ - Xdpp_) / Xd2_;
            Xd5_ = (Xdpp_ - Xl_) / Xd2_;
            Xq1_ = Xq_ - Xqp_;
            Xq2_ = Xqp_ - Xl_;
            Xq3_ = (Xqp_ - Xqpp_) / (Xq2_ * Xq2_);
            Xq4_ = (Xqp_ - Xqpp_) / Xq2_;
            Xq5_ = (Xqpp_ - Xl_) / Xq2_;
            Xqd_ = (Xq_ - Xl_) / (Xd_ - Xl_);
            G_ = Ra_ / (Ra_*Ra_ + Xqpp_*Xqpp_);
            B_ = -Xqpp_ / (Ra_*Ra_ + Xqpp_*Xqpp_);

        }

        ~GENROU() { }

        int allocate() override 
        { 
            f_.resize(size_);
            y_.resize(size_);
            yp_.resize(size_);
            tag_.resize(size_);
            fB_.resize(size_);
            yB_.resize(size_);
            ypB_.resize(size_);
            return 0; 
        }

        int initialize() override 
        { 
            double delta, omega, Eqp, psidp, psiqp, Edp, psiqpp, psidpp, psipp, 
                vd, vq, id, iq, ir, ii, Te, g, b, ksat;

            /* Initialization tricks -- assuming NO saturation */
            double vr, vi, p, q, vm2, Er, Ei;
            vr = Vr();
            vi = Vi();
            p = p0_;
            q = q0_;
            vm2 = vr*vr + vi*vi;
            Er = vr + (Ra_*p*vr + Ra_*q*vi - Xq_*p*vi + Xq_*q*vr) / vm2;
            Ei = vi + (Ra_*p*vi - Ra_*q*vr + Xq_*p*vr + Xq_*q*vi) / vm2;
            delta = atan2(Ei, Er);
            omega = 0;
            ir = (p*vr + q*vi) / vm2;
            ii = (p*vi - q*vr) / vm2;
            id = ir*sin(delta) -ii*cos(delta);
            iq = ir*cos(delta) +ii*sin(delta);
            vd = vr*sin(delta) - vi*cos(delta) + id*Ra_ - iq*Xqpp_;
            vq = vr*cos(delta) + vi*sin(delta) + id*Xqpp_ - iq*Ra_;
            psiqpp = -vd / (1 + omega);
            psidpp = vq / (1 + omega);
            Te = (psidpp - id*Xdpp_)*iq - (psiqpp - iq*Xdpp_)*id;
            psiqp = -(-(Xqp_-Xl_)*iq+psiqpp*(Xqp_-Xl_)/(Xqpp_-Xl_)) 
                /(1+(Xqp_-Xqpp_)/(Xqpp_-Xl_));
            Edp = psiqp - (Xqp_ - Xl_) * iq;
            psidp = -((Xdp_-Xl_)*id-psidpp*(Xdp_-Xl_)/(Xdpp_-Xl_))
                /(1+(Xdp_-Xdpp_)/(Xdpp_-Xl_));
            Eqp = psidp + (Xdp_ - Xl_) * id;

            /* Now we have the state variables, solve for alg. variables */

            y_[0] = delta; //= 0.55399038;
            y_[1] = omega; // = 0;
            y_[2] = Eqp; // = 0.995472581;
            y_[3] = psidp; // = 0.971299567;
            y_[4] = psiqp; // = 0.306880069;
            y_[5] = Edp; // = 0;

            y_[6] = psiqpp = -psiqp*Xq4_ - Edp*Xq5_;
            y_[7] = psidpp = psidp*Xd4_ + Eqp*Xd5_;
            y_[8] = psipp = sqrt(psiqpp*psiqpp + psidpp*psidpp);
            y_[9] = ksat = SB_*pow(psipp - SA_, 2);
            y_[10] = vd = -psiqpp*(1 + omega);
            y_[11] = vq = psidpp*(1 + omega);
            y_[12] = Te = (psidpp - id*Xdpp_)*iq - (psiqpp - iq*Xdpp_)*id; 
            y_[13] = id;
            y_[14] = iq;
            y_[15] = ir;
            y_[16] = ii;
            y_[17] = pmech_set_ = Te;
            y_[18] = efd_set_ = Eqp + Xd1_*(id 
                + Xd3_*(Eqp - psidp - Xd2_*id)) + psidpp*ksat;
            y_[19] = G_*(vd*sin(delta)+vq*cos(delta)) 
                - B_*(vd*-cos(delta) + vq*sin(delta)); /* inort, real */;
            y_[20] = B_*(vd*sin(delta)+vq*cos(delta)) 
                + G_*(vd*-cos(delta) + vq*sin(delta)); /* inort, imag */
            
            for (int i = 0; i < size_; ++i) yp_[i] = 0.0;
            return 0; 
        }

        int tagDifferentiable() override 
        { 
            for (int i = 0; i < size_; ++i) 
            {
                tag_[i] = i < 6;
            }
            return 0;
        }

        int evaluateResidual() override 
        { 
            double delta, omega, Eqp, psidp, psiqp, Edp, psiqpp, psidpp, psipp, 
                ksat, vd, vq, telec, id, iq, ir, ii, pmech, efd, inr, ini, vr, 
                vi, vr_inf, vi_inf, delta_dot, omega_dot, Eqp_dot, psidp_dot, 
                psiqp_dot, Edp_dot;

            /* Read variables */
            delta = y_[0];
            omega = y_[1];
            Eqp = y_[2];
            psidp = y_[3];
            psiqp = y_[4];
            Edp = y_[5];
            psiqpp = y_[6]; 
            psidpp = y_[7]; 
            psipp = y_[8]; 
            ksat = y_[9];
            vd = y_[10];
            vq = y_[11];
            telec = y_[12];
            id = y_[13];
            iq = y_[14];
            ir = y_[15];
            ii = y_[16];
            pmech = y_[17];
            efd = y_[18];
            inr = y_[19];
            ini = y_[20];
            vr = Vr(); 
            vi = Vi();
            
            /* Read derivatives */
            delta_dot = yp_[0];
            omega_dot = yp_[1];
            Eqp_dot = yp_[2];
            psidp_dot = yp_[3];
            psiqp_dot = yp_[4];
            Edp_dot = yp_[5];
    
            /* 6 GENROU differential equations */
            f_[0] = delta_dot - omega*(2*M_PI*60);
            f_[1] = omega_dot - (1/(2*H_)) * ((pmech - D_*omega) / (1 + omega) 
                - telec);
            f_[2] = Eqp_dot - (1/Tdop_) * (efd - (Eqp + Xd1_*(id + Xd3_*(Eqp 
                - psidp - Xd2_*id)) + psidpp*ksat ));
            f_[3] = psidp_dot - (1/Tdopp_) * (Eqp - psidp - Xd2_*id);
            f_[4] = psiqp_dot - (1/Tqopp_) * (Edp - psiqp + Xq2_*iq);
            f_[5] = Edp_dot - (1/Tqop_) * (-Edp + Xqd_*psiqpp*ksat 
                + Xq1_*(iq - Xq3_*(Edp + iq*Xq2_ - psiqp)));
            
            /* 11 GENROU algebraic equations */
            f_[6] = psiqpp - (-psiqp * Xq4_ - Edp * Xq5_);
            f_[7] = psidpp - (psidp * Xd4_ + Eqp * Xd5_);
            f_[8] = psipp - sqrt(pow(psidpp, 2.0) + pow(psiqpp, 2.0));
            f_[9] = ksat - SB_*pow(psipp - SA_, 2.0);
            f_[10] = vd + psiqpp * (1 + omega);
            f_[11] = vq - psidpp * (1 + omega);
            f_[12] = telec - ((psidpp - id*Xdpp_)*iq - (psiqpp - iq*Xdpp_)*id);
            f_[13] = id - (ir*sin(delta) - ii*cos(delta));
            f_[14] = iq - (ir*cos(delta) + ii*sin(delta));
            f_[15] = ir + G_*vr - B_*vi - inr;
            f_[16] = ii + B_*vr + G_*vi - ini;
            
            /* 2 GENROU control inputs are set to constant for this example */
            f_[17] = pmech - pmech_set_;
            f_[18] = efd - efd_set_;

            /* 2 GENROU current source definitions */
            f_[19] = inr - (G_*(sin(delta)*vd + cos(delta)*vq) 
                - B_*(-cos(delta)*vd + sin(delta)*vq));
            f_[20] = ini - (B_*(sin(delta)*vd + cos(delta)*vq) 
                + G_*(-cos(delta)*vd + sin(delta)*vq));

            /* Current balance */
            Ir() += inr - Vr()*G_ + Vi()*B_;
            Ii() += ini - Vr()*B_ - Vi()*G_;

            //printf("GENROU residual\n");
            //for (int i = 0 ; i < 21; ++i) printf("%d: %g\n", i, f_[i]);

            //printf("GENROU inr %g Vr %g B %g Vi %g G %g\n", inr, Vr(), B_, Vi(), G_);
            //printf("GENROU Ii = %g\n", inr - Vr()*B_ - Vi()*G_);

            return 0; 
        }

        int evaluateJacobian() override 
        { 
            /* TODO */
            return 0; 
        }

        /* Don't know what to do with any of these */
        int evaluateIntegrand() override { return 0; }
        int initializeAdjoint() override { return 0; }
        int evaluateAdjointResidual() override { return 0; }
        int evaluateAdjointIntegrand() override { return 0; }
        void updateTime(double t, double a) override {   }

    private:
        double& Vr()
        {
            return bus_->Vr();
        }

        double& Vi()
        {
            return bus_->Vi();
        }

        double& Ir()
        {
            return bus_->Ir();
        }

        double& Ii()
        {
            return bus_->Ii();
        }

    private:
    
        /* Identification */
        BaseBusT* bus_;
        const int busID_;
        int unit_id_;

        /* Initial terminal conditions */
        double p0_, q0_;

        /* Input parameters */
        double H_, D_, Ra_, Tdop_, Tdopp_, Tqopp_, Tqop_, Xd_, Xdp_, Xdpp_,
            Xq_, Xqp_, Xqpp_, Xl_, S10_, S12_;

        /* Derivied parameters */
        double SA_, SB_, Xd1_, Xd2_, Xd3_, Xd4_, Xd5_, Xq1_, Xq2_, Xq3_, Xq4_,
            Xq5_, Xqd_, G_, B_;

        /* Setpoints for control variables (determined at initialization) */
        double pmech_set_, efd_set_;
    };

}
}
