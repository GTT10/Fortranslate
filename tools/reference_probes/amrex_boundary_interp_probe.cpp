// Calls the pinned AMReX interpolation kernels; no copied kernel arithmetic.
#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_MFInterpolater.H>
#include <algorithm>
#include <array>
#include <iomanip>
#include <iostream>
#include <string>
int main(int argc, char** argv) {
    const bool wrong_time=argc>1 && std::string(argv[1])=="wrong_time";
    amrex::Initialize(argc, argv, false);
    {
        int cases; std::cin >> cases;
        std::cout << std::setprecision(17);
        for (int test=0; test<cases; ++test) {
            int id,axis,side,layout,profile,factor; double blend;
            std::cin >> id >> axis >> side >> layout >> profile >> factor >> blend;
            --axis;
            amrex::IntVect origin,shape,target;
            for(int d=0;d<3;++d) std::cin >> origin[d];
            for(int d=0;d<3;++d) std::cin >> shape[d];
            for(int d=0;d<3;++d) std::cin >> target[d];
            double root[4][4][4][8],future[4][4][4][8],lower[8],upper[8],actual[8];
            for(int k=0;k<4;++k) for(int j=0;j<4;++j) for(int i=0;i<4;++i)
                for(int c=0;c<8;++c) std::cin >> root[k][j][i][c];
            for(int k=0;k<4;++k) for(int j=0;j<4;++j) for(int i=0;i<4;++i)
                for(int c=0;c<8;++c) std::cin >> future[k][j][i][c];
            for(double& v:lower) std::cin >> v;
            for(double& v:upper) std::cin >> v;
            for(double& v:actual) std::cin >> v;
            AMREX_ALWAYS_ASSERT(std::cin.good());
            amrex::Box domain(amrex::IntVect(0),amrex::IntVect(3));
            amrex::Box fine_domain(amrex::IntVect(0),amrex::IntVect(4*factor-1));
            amrex::Box destination(origin,origin+shape-amrex::IntVect(1));
            destination.grow(1); destination &= fine_domain;
            amrex::IntVect ratio(factor);
            auto buffer=amrex::mf_lincc_interp.CoarseBox(destination,ratio);
            amrex::FArrayBox coarse(buffer,9),slope(amrex::coarsen(destination,ratio),27);
            amrex::FArrayBox result(amrex::Box(target,target),9);
            auto u=coarse.array();
            auto lower_index=buffer.smallEnd(),upper_index=buffer.bigEnd();
            for(int k=lower_index[2];k<=upper_index[2];++k)
                for(int j=lower_index[1];j<=upper_index[1];++j)
                    for(int i=lower_index[0];i<=upper_index[0];++i) {
                        std::array<int,3> point{i,j,k};
                        const double* old=root[std::clamp(k,0,3)][std::clamp(j,0,3)][std::clamp(i,0,3)];
                        const double* next=future[std::clamp(k,0,3)][std::clamp(j,0,3)][std::clamp(i,0,3)];
                        double fraction=blend;
                        if(point[axis]<0) { old=lower; next=lower; fraction=0; }
                        else if(point[axis]>3) { old=upper; next=upper; fraction=0; }
                        auto rhoe=[](const double* v) {
                            double vx=v[1]/v[0],vy=v[2]/v[0],vz=v[3]/v[0];
                            return v[4]-0.5*v[0]*(vx*vx+vy*vy+vz*vz);
                        };
                        for(int c=0;c<8;++c) u(i,j,k,c)=(1-fraction)*old[c]+fraction*next[c];
                        u(i,j,k,8)=(1-fraction)*rhoe(old)+fraction*rhoe(next);
                        // Deliberately wrong control: derive rhoe after blending.
                        if(wrong_time) {
                            double blended[8];
                            for(int c=0;c<8;++c) blended[c]=u(i,j,k,c);
                            u(i,j,k,8)=rhoe(blended);
                        }
                    }
            std::array<amrex::BCRec,9> bcs;
            for(auto& bc:bcs) for(int d=0;d<3;++d) {
                bc.setLo(d,d==axis?amrex::BCType::ext_dir:amrex::BCType::foextrap);
                bc.setHi(d,d==axis?amrex::BCType::ext_dir:amrex::BCType::foextrap);
            }
            auto donor=amrex::coarsen(target,ratio);
            amrex::mf_cell_cons_lin_interp_llslope(donor[0],donor[1],donor[2],slope.array(),
                coarse.const_array(),0,9,domain,ratio,bcs.data());
            for(int c=0;c<8;++c)
                amrex::mf_cell_cons_lin_interp(target[0],target[1],target[2],c,result.array(),0,
                    slope.const_array(),coarse.const_array(),0,9,ratio);
            std::cout << "CASE " << id;
            for(int c=0;c<8;++c) std::cout << ' ' << actual[c] << ' ' << result.const_array()(target[0],target[1],target[2],c);
            std::cout << '\n';
        }
    }
    amrex::Finalize();
}
