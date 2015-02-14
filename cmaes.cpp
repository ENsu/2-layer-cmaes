#include "cmaes.h"
#include "Eigen/Dense"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp" 


namespace Eigen {  
namespace internal {
template<typename Scalar>

struct scalar_normal_dist_op 
{
static boost::mt19937 rng;    // The uniform pseudo-random algorithm
mutable boost::normal_distribution<Scalar> norm;  // The gaussian combinator

EIGEN_EMPTY_STRUCT_CTOR(scalar_normal_dist_op)

    template<typename Index>
    inline const Scalar operator() (Index, Index = 0) const { return norm(rng); }
};

template<typename Scalar> boost::mt19937 scalar_normal_dist_op<Scalar>::rng;
template<typename Scalar>
struct functor_traits<scalar_normal_dist_op<Scalar> >
{ enum { Cost = 50 * NumTraits<Scalar>::MulCost, PacketAccess = false, IsRepeatable = false }; };

} // end namespace internal
} // end namespace Eigen
//from... http://stackoverflow.com/questions/6142576/sample-from-multivariate-normal-gaussian-distribution-in-c
