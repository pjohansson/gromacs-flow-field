#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/txtdump.h"

#include "gmx_io.h"

namespace flow
{

void read_rnemd_opts(std::vector<t_inpfile> &inp,
                     RNEMDOptions           &opts,
                     char                   *groups,
                     warninp                *wi)
{
    // Replace old "shear coupling" terminology with RNEMD
    replace_inp_entry(inp, "shear-coupling", "rnemd");
    replace_inp_entry(inp, "shear-axis", "rnemd-area-def-axis");
    replace_inp_entry(inp, "shear-direction", "rnemd-exchange-axis");
    replace_inp_entry(inp, "shear-strategy", "rnemd-strategy");
    replace_inp_entry(inp, "shear-tcoupl", "rnemd-tcoupl");
    replace_inp_entry(inp, "shear-area-size", "rnemd-area-size");
    replace_inp_entry(inp, "shear-zadj", "rnemd-zadj");
    replace_inp_entry(inp, "shear-ref-velocity", "rnemd-ref-velocity");
    replace_inp_entry(inp, "shear-grps", "rnemd-grps");

    printStringNewline(&inp, "REVERSE NON-EQUILIBRIUM MOLECULAR DYNAMICS (UNOFFICIAL)");

    printStringNoNewline(&inp, "Do RNEMD exchange of kinetic energy between areas");
    opts.bDoExchange = (getEnum<Boolean>(&inp, "rnemd", wi) != Boolean::No);

    printStringNoNewline(&inp, "Axis along which to create areas for the energy exchange ");
    printStringNoNewline(&inp, "and along which velocity vector to exchange energies: x, y or z");
    opts.area_def_axis = getEnum<flow::RnemdAreaDefAxis>(&inp, "rnemd-area-def-axis", wi);
    opts.energy_exchange_axis = getEnum<flow::RnemdEnergyExchangeAxis>(&inp, "rnemd-exchange-axis", wi);

    printStringNoNewline(&inp, "Strategy for setting up exchange areas: Edges or Edge-Center");
    printStringNoNewline(&inp, "Edges: exchange area 0 and 1 are respectively at the bottom and top");
    printStringNoNewline(&inp, "  edges of the system, along the selected axis");
    printStringNoNewline(&inp, "Edge-Center: exchange area 0 is split into the bottom and top edges");
    printStringNoNewline(&inp, "  of the system, area 1 is at the center");
    opts.strategy = getEnum<flow::RnemdStrategy>(&inp, "rnemd-strategy", wi);

    printStringNoNewline(&inp, "How often to perform the coupling");
    opts.tau    = get_ereal(&inp, "rnemd-tcoupl", 0.0, wi);

    printStringNoNewline(&inp, "Size of exchange areas and adjustment from the edges");
    opts.area_size = get_ereal(&inp, "rnemd-area-size", 0.0, wi);
    opts.zadj      = get_ereal(&inp, "rnemd-zadj", 0.0, wi);

    printStringNoNewline(&inp, "Reference velocity: Targeted velocity for both areas");
    printStringNoNewline(&inp, "  Area 0: -rnemd-ref-velocity");
    printStringNoNewline(&inp, "  Area 1: +rnemd-ref-velocity");
    opts.ref_velocity = get_ereal(&inp, "rnemd-ref-velocity", 0.0, wi);

    printStringNoNewline(&inp, "Groups to exchange for: must be 1 or 2, in the latter case ");
    printStringNoNewline(&inp, "for area 0 and 1 respectively");
    printStringNoNewline(&inp, "Note: this replaces user2-grps");
    setStringEntry(&inp, "rnemd-grps", groups, nullptr);
}


void check_rnemd_opts(const t_inputrec *ir,
                      warninp          *wi,
                      const gmx::EnumerationArray<PbcType, std::string> pbcTypeNames)
{
    const auto& opts = ir->rnemd_opts;

    if (opts.bDoExchange)
    {
        if ((opts.strategy == flow::RnemdStrategy::Edges)
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
    serializer->doEnumAsInt(&opts.area_def_axis);
    serializer->doEnumAsInt(&opts.energy_exchange_axis);
    serializer->doEnumAsInt(&opts.strategy);
    serializer->doReal(&opts.tau);
    serializer->doReal(&opts.area_size);
    serializer->doReal(&opts.zadj);
    serializer->doReal(&opts.ref_velocity);
}

void pr_rnemd(FILE* fp, int indent, const RNEMDOptions &opts)
{
    pr_str(fp, indent, "rnemd", booleanValueToString(opts.bDoExchange));
    pr_str(fp, indent, "rnemd-area-def-axis", enumValueToString(opts.area_def_axis));
    pr_str(fp, indent, "rnemd-exchange-axis", enumValueToString(opts.energy_exchange_axis));
    pr_str(fp, indent, "rnemd-strategy", enumValueToString(opts.strategy));
    pr_real(fp, indent, "rnemd-tcoupl", opts.tau);
    pr_real(fp, indent, "rnemd-area-size", opts.area_size);
    pr_real(fp, indent, "rnemd-zadj", opts.zadj);
    pr_real(fp, indent, "rnemd-ref-velocity", opts.ref_velocity);
}

} // namespace flow
