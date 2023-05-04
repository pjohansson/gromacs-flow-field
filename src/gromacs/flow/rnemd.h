#ifndef FLOW_RNEMD_H
#define FLOW_RNEMD_H

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

using namespace flow;

struct RNEMD {
    FILE *log_pexchange;         /* Opened file to log momentum exchange into */

    RnemdAreaDefAxis area_def_axis;                 /* Axis along which the areas that are coupled are defined */
    RnemdEnergyExchangeAxis energy_exchange_axis;   /* Axis for which the kinetic energy itself is swapped along */

    RnemdStrategy strategy;  /* Strategy for determining the exchange areas */

    size_t num_groups,      /* Number of groups to include atoms from
                               Will be 1 or 2, if 1 both areas will include atoms from
                               the group, but if 2, area 0 and 1 will include atoms only
                               from groups 0 and 1, respectively                      */

           step;            /* Couple the velocities at multiples of this step */

    real   size,            /* Area size along the axis */
           zedge_adj,       /* Adjust the zmin position by this amount from the edges */
           ref_velocity;    /* Reference velocity to target in the areas              */
                            /* Area 0: velocity = -ref_velocity                       */
                            /* Area 1: velocity = +ref_velocity                       */
};

RNEMD init_rnemd(const t_inputrec              *ir,
                 const matrix                   box,
                 const SimulationGroups        *groups,
                 const t_commrec               *cr,
                 const char                    *fnlog,
                 const struct gmx_output_env_t *oenv,
                 const gmx::MDLogger           &mdlog);

void do_rnemd_exchange(t_state                *state,
                       const t_mdatoms        *mdatoms,
                       const int64_t           current_step,
                       const double            current_time,
                       const RNEMD            &rnemd,
                       const SimulationGroups *groups,
                       const t_commrec        *cr);

#endif // FLOW_RNEMD_H