# Public problem names and their source directories.
#
# This manifest is intentionally free of build rules.  It is the single source
# of truth that can be included by both the repository root and
# problems/GNUmakefile without duplicating the public-name/directory mapping.
PROBLEMS := l2 heat eriksson pure_diffusion_igrm oil igrm_l2 igrm_heat \
	igrm_eirksson igrm_stokes igrm_pollution

l2_DIR := L2Projection
heat_DIR := heat
eriksson_DIR := eriksson
pure_diffusion_igrm_DIR := pure_diffusion_igrm
oil_DIR := oil
igrm_l2_DIR := igrm_L2Projection
igrm_heat_DIR := igrm_heat
igrm_eirksson_DIR := igrm_eirksson
igrm_stokes_DIR := igrm_stokes
igrm_pollution_DIR := igrm_pollution
