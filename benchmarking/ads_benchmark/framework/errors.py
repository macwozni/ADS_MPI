"""Framework-specific exceptions with concise command-line diagnostics."""


class BenchmarkError(Exception):
    """Base class for expected benchmark configuration/runtime failures."""


class ConfigurationError(BenchmarkError):
    """A JSON profile is malformed or incomplete."""


class ValidationError(BenchmarkError):
    """An expanded benchmark case violates a public invariant."""


class RegistryError(BenchmarkError):
    """A component name is invalid, unknown, or registered twice."""


class DuplicateCaseError(BenchmarkError):
    """A profile expands to the same semantic case more than once."""


class ProvenanceError(BenchmarkError):
    """Repository provenance cannot be established reliably."""


class StorageError(BenchmarkError):
    """A result path is unsafe, already exists, or cannot be written."""


class ExecutionError(BenchmarkError):
    """A registered component cannot execute the requested case."""
