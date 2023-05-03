#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/txtdump.h"

#include "gmx_io.h"

namespace flow
{

void read_rnemd_opts(std::vector<t_inpfile> &inp,
                     RNEMDOptions           &opts,
                     char                   *groups,
                     warninp                *wi)
{
    printStringNewline(&inp, "REVERSE NON-EQUILIBRIUM MOLECULAR DYNAMICS (UNOFFICIAL)");

    printStringNoNewline(&inp, "Do RNEMD exchange of kinetic energy between areas");
    opts.bDoExchange = (getEnum<Boolean>(&inp, "shear-coupling", wi) != Boolean::No);

    printStringNoNewline(&inp, "Axis along which to exchange the energies and the direction");
    printStringNoNewline(&inp, "along which to shear: x, y or z");
    opts.axis = getEnum<flow::ShearAxis_axis>(&inp, "shear-axis", wi);
    opts.direction = getEnum<flow::ShearAxis_direction>(&inp, "shear-direction", wi);

    printStringNoNewline(&inp, "Strategy for setting up exchange areas: Edges or Edge-Center");
    printStringNoNewline(&inp, "Edges: exchange area 0 and 1 are respectively at the bottom and top");
    printStringNoNewline(&inp, "  edges of the system, along the selected axis");
    printStringNoNewline(&inp, "Edge-Center: exchange area 0 is split into the bottom and top edges");
    printStringNoNewline(&inp, "  of the system, area 1 is at the center");
    opts.strategy = getEnum<flow::ShearCouplStrategy>(&inp, "shear-strategy", wi);

    printStringNoNewline(&inp, "How often to perform the coupling");
    opts.tau    = get_ereal(&inp, "shear-tcoupl", 0.0, wi);

    printStringNoNewline(&inp, "Size of exchange areas and adjustment from the edges");
    opts.area_size = get_ereal(&inp, "shear-area-size", 0.0, wi);
    opts.zadj      = get_ereal(&inp, "shear-zadj", 0.0, wi);

    printStringNoNewline(&inp, "Reference velocity: Targeted velocity for both areas");
    printStringNoNewline(&inp, "  Area 0: -shear-ref-velocity");
    printStringNoNewline(&inp, "  Area 1: +shear-ref-velocity");
    opts.ref_velocity = get_ereal(&inp, "shear-ref-velocity", 0.0, wi);

    printStringNoNewline(&inp, "Groups to shear with: must be 1 or 2, in the latter case ");
    printStringNoNewline(&inp, "for area 0 and 1 respectively");
    printStringNoNewline(&inp, "Note: this replaces user2-grps");
    setStringEntry(&inp, "shear-grps", groups, nullptr);
}


void check_rnemd_opts(const t_inputrec *ir,
                      warninp          *wi,
                      const gmx::EnumerationArray<PbcType, std::string> pbcTypeNames)
{
    const auto& opts = ir->rnemd_opts;

    if (opts.bDoExchange)
    {
        if ((opts.strategy == flow::ShearCouplStrategy::Edges)
            && (ir->pbcType != PbcType::XY))
        {
            char warn_buf[STRLEN];

            sprintf(warn_buf,
                "With shear-strategy = %s the system should probably not "
                "be periodic along the axis, since both edges will shear "
                "against each other. Use pbc = %s unless you are sure. "
                "(currently pbc = %s)",
                enumValueToString(opts.strategy),
                pbcTypeNames[PbcType::XY].c_str(),
                pbcTypeNames[ir->pbcType].c_str()
            );
            warning(wi, warn_buf);
        }
    }
}

//! Verify that flow field options are correctly set, after all groups have been set
void check_rnemd_groups(const t_inputrec *ir,
                        const std::vector<std::string> &group_names,
                        warninp   *wi)
{
    const auto& opts = ir->rnemd_opts;

    if (opts.bDoExchange
        && (group_names.empty() || group_names.size() > 2))
    {
        char warn_buf[STRLEN];

        snprintf(warn_buf, STRLEN,
            "Invalid shear-grps input: must be 1 or 2 groups (is %lu), "
            "in which case area 0 only couples to atoms in the first "
            "group and area 1 to atoms in the second.",
            group_names.size()
        );

        warning_error(wi, warn_buf);
    }
}

void do_tpx_rnemd(gmx::ISerializer *serializer,
                  RNEMDOptions     &opts)
{
    serializer->doBool(&opts.bDoExchange);
    serializer->doEnumAsInt(&opts.axis);
    serializer->doEnumAsInt(&opts.direction);
    serializer->doEnumAsInt(&opts.strategy);
    serializer->doReal(&opts.tau);
    serializer->doReal(&opts.area_size);
    serializer->doReal(&opts.zadj);
    serializer->doReal(&opts.ref_velocity);
}


void pr_rnemd(FILE* fp, int indent, const RNEMDOptions &opts)
{
    pr_str(fp, indent, "shear-coupling", booleanValueToString(opts.bDoExchange));
    pr_str(fp, indent, "shear-axis", enumValueToString(opts.axis));
    pr_str(fp, indent, "shear-direction", enumValueToString(opts.direction));
    pr_str(fp, indent, "shear-strategy", enumValueToString(opts.strategy));
    pr_real(fp, indent, "shear-tcoupl", opts.tau);
    pr_real(fp, indent, "shear-area-size", opts.area_size);
    pr_real(fp, indent, "shear-zadj", opts.zadj);
    pr_real(fp, indent, "shear-ref-velocity", opts.ref_velocity);
}

} // namespace flow
