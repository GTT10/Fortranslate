// Serial indexing/storage shim for unchanged pinned CTU bodies.
// Unsupported scenario branches terminate instead of approximating them.
#include <array>
#include <vector>
#include <type_traits>
#include <stdexcept>
#define AMREX_D_DECL(a,b,c) a,b,c
#define AMREX_D_TERM(a,b,c) a b c
#define AMREX_D_PICK(a,b,c) c
#define AMREX_GPU_DEVICE
#define AMREX_SPACEDIM 3
#define NUM_ADV 0
#define NUM_LIN 0
constexpr int QRHO=0,QU=1,QV=2,QW=3,QPRES=4,QFS=5,QREINT=12,QVAR=13;
constexpr int URHO=0,UMX=1,UMY=2,UMZ=3,UEDEN=4,UFS=5,UEINT=12,UTEMP=13,NVAR=14;
constexpr int GDU=0,GDV=1,GDW=2,GDPRES=3,GDGAME=4,NGDNV=5,QC=0,QGAMC=1;
constexpr int im2=0,im1=1,i0=2,ip1=3,ip2=4;
namespace PCPhysBCType { enum { interior, inflow, user_bc, outflow, slip_wall, no_slip_wall, symmetry }; }
namespace amrex {
namespace Math { using std::copysign; }
namespace PhysBCType { constexpr int noslipwall=PCPhysBCType::no_slip_wall; }
template<class T,int N> using GpuArray=std::array<T,N>;
struct IntVect {
 int v[3];
 int& operator[](int n){return v[n];} int operator[](int n)const{return v[n];}
 static IntVect TheDimensionVector(int n){return {{n==0,n==1,n==2}};}
 IntVect operator+(IntVect b)const{return {{v[0]+b[0],v[1]+b[1],v[2]+b[2]}};}
 IntVect operator-(IntVect b)const{return {{v[0]-b[0],v[1]-b[1],v[2]-b[2]}};}
};
IntVect operator*(int a,IntVect b){return {{a*b[0],a*b[1],a*b[2]}};}
struct Box {
 IntVect lo,hi;bool node[3]={false,false,false};
 Box& grow(int d,int n){lo[d]-=n;hi[d]+=n;return *this;}
 IntVect loVect3d()const{return lo;} IntVect hiVect3d()const{return hi;}
 void setSmall(int d,int n){lo[d]=n;} void setBig(int d,int n){hi[d]=n;}
 bool ok()const{return lo[0]<=hi[0]&&lo[1]<=hi[1]&&lo[2]<=hi[2];}
};
Box grow(Box b,int n){for(int d=0;d<3;d++){b.lo[d]-=n;b.hi[d]+=n;}return b;}
Box grow(Box b,int d,int n){b.lo[d]-=n;b.hi[d]+=n;return b;}
Box growHi(Box b,int d,int n){b.hi[d]+=n;return b;}
Box surroundingNodes(Box b,int d){b.hi[d]++;b.node[d]=true;return b;}
Box enclosedCells(Box b){for(int d=0;d<3;d++)if(b.node[d]){b.hi[d]--;b.node[d]=false;}return b;}
template<class F>void ParallelFor(Box b,F f){for(int k=b.lo[2];k<=b.hi[2];k++)for(int j=b.lo[1];j<=b.hi[1];j++)for(int i=b.lo[0];i<=b.hi[0];i++)f(i,j,k);}
void* The_Async_Arena(){return nullptr;}
[[noreturn]]void Error(const char* s){throw std::runtime_error(s);}
template<class T>struct Array4 {
 T* p=nullptr; Box b{};int nc=0;
 Array4()=default;Array4(T* ptr,Box box,int n):p(ptr),b(box),nc(n){}
 template<class U,typename=std::enable_if_t<std::is_convertible_v<U*,T*>>>
 Array4(const Array4<U>& a):p(a.p),b(a.b),nc(a.nc){}
 explicit operator bool()const{return p!=nullptr;}
 bool contains(int i,int j,int k)const{return i>=b.lo[0]&&i<=b.hi[0]&&j>=b.lo[1]&&j<=b.hi[1]&&k>=b.lo[2]&&k<=b.hi[2];}
 T& operator()(int i,int j,int k,int n=0)const{
  if(!p||i<b.lo[0]||i>b.hi[0]||j<b.lo[1]||j>b.hi[1]||k<b.lo[2]||k>b.hi[2]||n<0||n>=nc)Error("array bounds");
  return p[(((k-b.lo[2])*(b.hi[1]-b.lo[1]+1)+j-b.lo[1])*(b.hi[0]-b.lo[0]+1)+i-b.lo[0])*nc+n];
 }
 T& operator()(IntVect v,int n=0)const{return (*this)(v[0],v[1],v[2],n);}
};
struct FArrayBox {
 Box b;int nc;std::vector<double> data;
 FArrayBox(Box box,int n,void* =nullptr):b(box),nc(n),data((b.hi[0]-b.lo[0]+1)*(b.hi[1]-b.lo[1]+1)*(b.hi[2]-b.lo[2]+1)*n,0.){}
 Array4<double> array(){return {data.data(),b,nc};}
};
}
template<class... A>void weno_reconstruct_5js(A...){amrex::Error("WENO outside oracle scope");}
template<class... A>void weno_reconstruct_5z(A...){amrex::Error("WENO outside oracle scope");}
template<class... A>void weno_reconstruct_7z(A...){amrex::Error("WENO outside oracle scope");}
template<class... A>void weno_reconstruct_3z(A...){amrex::Error("WENO outside oracle scope");}
