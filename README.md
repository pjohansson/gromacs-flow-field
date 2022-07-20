# Gromacs Flow Field

This fork of [Gromacs](http://www.gromacs.org/) modifies it to enable output
of two-dimensional flow fields, collected during a simulation.

## Installation

Follow the regular Gromacs [installations instructions](http://manual.gromacs.org/current/install-guide/index.html).
While not required, it is recommended to set the `cmake` option
`-DGMX_VERSION_STRING_OF_FORK=flow-field`.

## Usage

The flow field collection uses the following [MDP options](http://manual.gromacs.org/documentation/current/user-guide/mdp-options.html):

```mdp
; Do flow field collection: No or Yes
flow-field               = yes
; This selects the subset of atoms for the flow field
; collection. You can select multiple groups, in which
; case the fields for all groups combined and the fields
; for all individual groups are all collected and written
; to disk. If no groups are selected, all atoms in the
; system will be used.
flow-field-grps          = water glycerol
; Interval in steps between sampling flow field data
flow-nstsample           = 10
; Interval in steps between averaging and outputting flow field data
flow-nstoutput           = 5000
; Number of flow field grid bins along x
flow-nx                  = 100
; Number of flow field grid bins along z
flow-nz                  = 110
```

Additionally, `mdrun` needs to know where to save the output, for which the `-flow` option has been added:

```bash
$ gmx mdrun -flow maps/flow   # Saves to `maps/flow_00001.dat`,
                              #          `maps/flow_00002.dat`,
                              #          `maps/flow_00003.dat`,
                              # etc.
```

## Limitations

*   Can currently only sample data in the x-z plane.
*   Temperature calculation is only correct for water
*   Only works for static box size, so no pressure coupling

## File formats

Tools used to read and manipulate the created data is available in a Python 3 module.
It is available in this repository:

[https://github.com/pjohansson/gmx_flow_utils](https://github.com/pjohansson/gmx_flow_utils)

## License

These changes are distributed under the same terms as Gromacs. See
[COPYING](COPYING) for more information.
