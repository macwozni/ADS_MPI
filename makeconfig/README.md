# Example make configurations

The repository-root `m_options` is the active local configuration. Files in
this directory are selectable examples; they reuse its dependency paths and
override the compiler/build-mode settings for a particular case.

```bash
make CONFIG=makeconfig/gnu-debug.mk show-config
make CONFIG=makeconfig/gnu-release.mk all
```

Edit the shared dependency paths in `m_options` and select an example directly
with `CONFIG=...`. If you create another overlay, keep it in this directory or
adjust its relative include of root `m_options`.

The supported scope is GCC/GFortran with MUMPS. Intel compiler overlays and
ParMETIS are not qualified paths for current work.
