#include <iostream>
#define CL_TARGET_OPENCL_VERSION 120
#include <boost/compute/core.hpp>
#include <boost/compute/algorithm/copy.hpp>
#include <boost/compute/algorithm/transform.hpp>
#include <boost/compute/container/vector.hpp>
#include <boost/compute/types/fundamental.hpp>

namespace compute = boost::compute;
using compute::float4_;

struct Particle {
    float4_ position;
    float4_ velocity;
    float density;
    float time;
};
BOOST_COMPUTE_ADAPT_STRUCT(Particle, Particle, (position, velocity, time))
BOOST_COMPUTE_FUNCTION(Particle, step, (Particle particle), {
    // determine possible time step dt
    // reconstruct neighbor attributes
    // compute density
    // compute forces
    // compute velocity using forces
    // compute new density
    // compute pressure and pressure forces
    // integrate particle over time using dt
});

int main()
{
    compute::device device = compute::system::default_device();
    compute::context context(device);
    compute::command_queue queue(context, device);


    compute::vector<float4_> input(5, context);
    compute::vector<float4_> output(5, context);
    compute::copy(
        reinterpret_cast<Particle *>(points),
        reinterpret_cast<Particle *>(points) + 5,
        input.begin(),
        queue
    );
    compute::transform(input.begin(), input.end(), output.begin(), step, queue);

    return 0;
}