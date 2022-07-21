#include "accelerate.h"

namespace flow
{

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
