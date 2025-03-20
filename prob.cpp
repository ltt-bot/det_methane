#include "prob.H"
#include "mechanism.H"
#include "IndexDefines.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  amrex::ParmParse pp("prob");

  // Need to fill all the values in ProbParmDevice
  pp.query("equiv_ratio", PeleC::h_prob_parm_device->equiv_ratio);
  pp.query("reac_temp", PeleC::h_prob_parm_device->reac_temp);
  pp.query("reac_pres", PeleC::h_prob_parm_device->reac_pres);

  pp.query("inj_p0", PeleC::h_prob_parm_device->inj_p0);
  pp.query("inj_t0", PeleC::h_prob_parm_device->inj_t0);
  pp.query("inject_fuel", PeleC::h_prob_parm_device->inject_fuel);
  pp.query("num_inlet", PeleC::h_prob_parm_device->num_inlet);
  pp.query("inj_R", PeleC::h_prob_parm_device->inj_R);
  pp.query("p_outlet", PeleC::h_prob_parm_device->p_outlet);
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}