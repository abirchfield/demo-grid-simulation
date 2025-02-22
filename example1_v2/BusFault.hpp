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

    class BusFault : public ComponentT
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
        BusFault(BaseBusT* bus) : bus_(bus), R_(0), X_(0.01), 
            status_(0), busID_(0)
        {
            size_ = 0;
        }

        BusFault(BaseBusT* bus, double R, double X, int status) : 
            bus_(bus), R_(R), X_(X), status_(status), busID_(0)
        {
            size_ = 0;
        }

        ~BusFault() { }

        int allocate() override { return 0; }

        int initialize() override { return 0; }

        int tagDifferentiable() override { return 0; }

        int evaluateResidual() override 
        { 
            if (status_)
            {
                double B = -X_ / (X_*X_ + R_*R_);
                double G = R_ / (X_*X_ + R_*R_);
                Ir() += -Vr()*G + Vi()*B;
                Ii() += -Vr()*B - Vi()*G;
            }
            return 0; 
        }

        int evaluateJacobian() override { return 0; }


        int evaluateIntegrand() override { return 0; }
        int initializeAdjoint() override { return 0; }
        int evaluateAdjointResidual() override { return 0; }
        int evaluateAdjointIntegrand() override { return 0; }

        void updateTime(double t, double a) override {   }

    public:
        void setR(double R)  { R_ = R; }

        void setX(double X)  { X_ = X; }

        void setStatus(int status) { status_ = status; }


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
        BaseBusT* bus_;
        double R_;
        double X_;
        int status_;
        const int busID_;
    };

}
}
