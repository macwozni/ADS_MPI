# Example make configurations

The repository-root `m_options` is the active local configuration. Files in
this directory are selectable examples; they reuse its dependency paths and
override the compiler/build-mode settings for a particular case.

```bash
make CONFIG=makeconfig/gnu-debug.mk show-config
make CONFIG=makeconfig/gnu-release.mk all
make CONFIG=makeconfig/intel-oneapi-debug.mk all
make CONFIG=makeconfig/intel-oneapi-release.mk all
```

Edit the shared dependency paths in `m_options` and select an example directly
with `CONFIG=...`. If you create another overlay, keep it in this directory or
adjust its relative include of root `m_options`.

The Intel examples cover library and problem builds. Running the complete test
tree with Intel additionally requires compiler-neutral test flags throughout
the individual suites and an Intel-built pFUnit; the current suite Makefiles
use GNU-specific diagnostic flags. All numerical libraries must also be
compatible with the selected Intel MPI compiler wrapper.
