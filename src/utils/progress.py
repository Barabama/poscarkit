import sys
from typing import Iterable, Iterator, Optional, Any, Generator
from contextlib import contextmanager

try:
    from tqdm import tqdm as TqdmProgress
except ImportError:
    TqdmProgress = None


class SimpleProgress:
    """
    A simple progress bar wrapper around tqdm.
    Supports both fixed-length and variable-length iterables.
    Display width is limited to 80 characters.
    """

    def __init__(
        self,
        iterable: Optional[Iterable] = None,
        total: Optional[int] = None,
        desc: str = "",
        unit: str = "it",
        unit_scale: bool = False,
        ncols: int = 80,
        mininterval: float = 0.1,
        dynamic_ncols: bool = False,
        leave: bool = True,
    ):
        """
        Initialize progress bar.

        Args:
            iterable: Iterable to wrap
            total: Total number of items (None for unknown length)
            desc: Description prefix
            unit: Unit string
            unit_scale: Whether to scale units (K, M, G)
            ncols: Display width in characters
            mininterval: Minimum update interval in seconds
            dynamic_ncols: Whether to dynamically adjust width
            leave: Whether to leave progress bar after completion
        """
        if TqdmProgress is None:
            raise ImportError("tqdm is required. Install it with: pip install tqdm")

        self.iterable = iterable
        self.total = total
        self.desc = desc
        self.unit = unit
        self.unit_scale = unit_scale
        self.ncols = ncols
        self.mininterval = mininterval
        self.dynamic_ncols = dynamic_ncols
        self.leave = leave

        self.n = 0
        self._tqdm: Optional[TqdmProgress] = None

        if iterable is not None:
            if total is None:
                try:
                    total = len(iterable)  # type: ignore
                except (TypeError, AttributeError):
                    total = None
            self.total = total

    def _create_tqdm(self) -> TqdmProgress:
        """Create tqdm instance."""
        return TqdmProgress(
            iterable=self.iterable,
            total=self.total,
            desc=self.desc,
            unit=self.unit,
            unit_scale=self.unit_scale,
            ncols=self.ncols,
            mininterval=self.mininterval,
            dynamic_ncols=self.dynamic_ncols,
            leave=self.leave,
            file=sys.stdout,
        )

    def update(self, n: int = 1):
        """Update progress by n steps."""
        if self._tqdm is None:
            self._tqdm = self._create_tqdm()
        self._tqdm.update(n)
        self.n += n

    def write(self, msg: str = "", end: str = "\n"):
        """
        Write a message to stdout, ensuring progress bar is on a new line.
        
        Args:
            msg: Message to write
            end: String to append at the end (default: newline)
            
        Examples:
            >>> with progress(total=100, desc="Processing") as p:
            ...     for i in range(100):
            ...         p.write(f"Processing item {i}")
            ...         p.update(1)
        """
        if self._tqdm is None:
            self._tqdm = self._create_tqdm()
        
        self._tqdm.write(msg, end=end)

    def wrap(self, iterable: Optional[Iterable] = None) -> Iterator[Any]:
        """
        Wrap an iterable to automatically update progress.

        Args:
            iterable: Iterable to wrap (uses self.iterable if not provided)

        Yields:
            Items from iterable

        Examples:
            >>> for item in SimpleProgress(desc="Processing").wrap(data):
            ...     process(item)
        """
        if iterable is None:
            iterable = self.iterable

        if iterable is None:
            raise ValueError("No iterable provided")

        if self._tqdm is None:
            if self.total is None:
                try:
                    self.total = len(iterable)  # type: ignore
                except (TypeError, AttributeError):
                    pass
            self._tqdm = self._create_tqdm()

        try:
            for item in iterable:
                self.n += 1
                self._tqdm.update(1)
                yield item
        finally:
            self.close()

    def close(self):
        """Close progress bar."""
        if self._tqdm is not None:
            self._tqdm.close()
            self._tqdm = None

    def __iter__(self) -> Iterator[Any]:
        """Iterate over wrapped iterable."""
        if self._tqdm is None:
            self._tqdm = self._create_tqdm()

        try:
            for item in self._tqdm:
                self.n += 1
                yield item
        finally:
            self.close()

    def __enter__(self):
        """Context manager entry."""
        if self._tqdm is None:
            self._tqdm = self._create_tqdm()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()
        return False


def progress(
    iterable: Optional[Iterable] = None,
    total: Optional[int] = None,
    desc: str = "",
    unit: str = "it",
    unit_scale: bool = False,
    ncols: int = 80,
    mininterval: float = 0.1,
    dynamic_ncols: bool = False,
    leave: bool = True,
) -> SimpleProgress:
    """
    Create a progress bar.

    Args:
        iterable: Iterable to wrap
        total: Total number of items (None for unknown length)
        desc: Description prefix
        unit: Unit string
        unit_scale: Whether to scale units (K, M, G)
        ncols: Display width in characters
        mininterval: Minimum update interval in seconds
        dynamic_ncols: Whether to dynamically adjust width
        leave: Whether to leave progress bar after completion

    Returns:
        SimpleProgress instance

    Examples:
        >>> for i in progress(range(100), desc="Processing"):
        ...     time.sleep(0.01)

        >>> with progress(total=None, desc="Yield") as p:
        ...     for item in generator():
        ...         p.write(f"Processing {item}")
        ...         p.update(1)
    """
    return SimpleProgress(
        iterable=iterable,
        total=total,
        desc=desc,
        unit=unit,
        unit_scale=unit_scale,
        ncols=ncols,
        mininterval=mininterval,
        dynamic_ncols=dynamic_ncols,
        leave=leave,
    )


@contextmanager
def progress_context(
    total: Optional[int] = None,
    desc: str = "",
    unit: str = "it",
    unit_scale: bool = False,
    ncols: int = 80,
    mininterval: float = 0.1,
    dynamic_ncols: bool = False,
    leave: bool = True,
) -> Generator[SimpleProgress, None, None]:
    """
    Context manager for progress bar with unknown length.

    Args:
        total: Total number of items (None for unknown length)
        desc: Description prefix
        unit: Unit string
        unit_scale: Whether to scale units (K, M, G)
        ncols: Display width in characters
        mininterval: Minimum update interval in seconds
        dynamic_ncols: Whether to dynamically adjust width
        leave: Whether to leave progress bar after completion

    Yields:
        SimpleProgress instance

    Examples:
        >>> with progress_context(desc="Processing") as p:
        ...     for item in generator():
        ...         p.write(f"Processing {item}")
        ...         p.update(1)
    """
    p = SimpleProgress(
        iterable=None,
        total=total,
        desc=desc,
        unit=unit,
        unit_scale=unit_scale,
        ncols=ncols,
        mininterval=mininterval,
        dynamic_ncols=dynamic_ncols,
        leave=leave,
    )
    try:
        yield p
    finally:
        p.close()
