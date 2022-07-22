#include "accelerate.h"

namespace flow
{

//! Calculate the current acceleration multiplier
//!
//! Uses a smooth-step function to slowly increase the acceleration
//! multiplier from 0.0 to 1.0. If step_complete <= 0, always returns 1.0.
real calc_acceleration_multiplier(const int64_t step,
                                  const int64_t step_complete)
{
    if ((step >= step_complete) || (step_complete <= 0))
    {
        return 1.0;
    }
    else
    {
        const auto x = static_cast<real>(step) / static_cast<real>(step_complete);

        // Smoothstep function for range [0.0, 1.0) -> [0.0, 1.0)
        return 6.0 * powf(x, 5.0) - 15.0 * powf(x, 4.0) + 10 * powf(x, 3.0);
    }
}

void print_local_acceleration_info(const flow::LocalAcceleration &opts,
                                   const gmx::MDLogger           &mdlog)
{
    // Log to warning level, which prints both to md.log and stdout
    // (info level only writes to md.log)

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("**********************************\n")
        .appendText("* LOCAL ACCELERATION INFORMATION *\n")
        .appendText("**********************************");

    const auto& rmin = opts.rmin;
    const auto& rmax = opts.rmax;
    const auto& axis = opts.check_axis;

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Local acceleration is enabled: %s\n",
            opts.doLocalAcceleration ? "yes" : "no"
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Origin:  [%g, %g, %g]\n",
            axis[XX] ? rmin[XX] : -1,
            axis[YY] ? rmin[YY] : -1,
            axis[ZZ] ? rmin[ZZ] : -1
        )
        .appendTextFormatted(
            "End:     [%g, %g, %g]\n",
            axis[XX] ? rmax[XX] : -1,
            axis[YY] ? rmax[YY] : -1,
            axis[ZZ] ? rmax[ZZ] : -1
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendTextFormatted(
            "Activation time: %g",
            opts.tau
        );

    GMX_LOG(mdlog.warning)
        .asParagraph()
        .appendText("**************************************\n")
        .appendText("* END LOCAL ACCELERATION INFORMATION *\n")
        .appendText("**************************************");
}

} // namespace flow
