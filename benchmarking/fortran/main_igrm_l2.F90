program igrm_l2_manufactured

   use benchmark_contract, ONLY: BenchmarkAdapter
   use benchmark_harness, ONLY: RunManufacturedBenchmark
   use igrm_l2_benchmark_adapter, ONLY: RegisterIgrmL2Adapter

   implicit none

   type(BenchmarkAdapter) :: adapter

   call RegisterIgrmL2Adapter(adapter)
   call RunManufacturedBenchmark(adapter)

end program igrm_l2_manufactured
