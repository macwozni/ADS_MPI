program igrm_heat_manufactured

   use benchmark_contract, ONLY: BenchmarkAdapter
   use benchmark_harness, ONLY: RunManufacturedBenchmark
   use igrm_heat_benchmark_adapter, ONLY: RegisterIgrmHeatAdapter

   implicit none

   type(BenchmarkAdapter) :: adapter

   call RegisterIgrmHeatAdapter(adapter)
   call RunManufacturedBenchmark(adapter)

end program igrm_heat_manufactured
