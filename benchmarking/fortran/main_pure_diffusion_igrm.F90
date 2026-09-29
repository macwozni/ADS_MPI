program pure_diffusion_igrm_manufactured

   use benchmark_contract, ONLY: BenchmarkAdapter
   use benchmark_harness, ONLY: RunManufacturedBenchmark
   use pure_diffusion_igrm_benchmark_adapter, ONLY: &
      RegisterPureDiffusionIgrmAdapter

   implicit none

   type(BenchmarkAdapter) :: adapter

   call RegisterPureDiffusionIgrmAdapter(adapter)
   call RunManufacturedBenchmark(adapter)

end program pure_diffusion_igrm_manufactured
