from typing import Tuple
import os

# IO module for read_breaker
class FastqWriter:
    """
    Writes paired-end reads to gzipped FASTQ files.

    Attributes:
        read1_path (str): Output path for R1 FASTQ file.
        read2_path (str): Output path for R2 FASTQ file.
    """

    def __init__(self, output_dir: str, filename_stub: str, mode: str = "w") -> None:
        """
        Initializes the FastqWriter instance.

        Args:
            output_dir (str): Directory where FASTQ files will be saved.
            filename_stub (str): Base name for output files.
            mode (str): File mode ('w' for write, 'a' for append).
        """
        mode += "t"  # ensure text mode for gzip
        self.read1_path = os.path.join(output_dir, f"{filename_stub}.R1.fastq.gz")
        self.read2_path = os.path.join(output_dir, f"{filename_stub}.R2.fastq.gz")
        self.read1_file = gzip.open(self.read1_path, mode, compresslevel=3)
        self.read2_file = gzip.open(self.read2_path, mode, compresslevel=3)

    def write(self, read_pair: Tuple[str, str, str, str, str]) -> None:
        """
        Writes a single paired-end read to the FASTQ files.

        Args:
            read_pair: Tuple of (read_id, seq1, qual1, seq2, qual2).
        """
        read_id, seq1, qual1, seq2, qual2 = read_pair
        self.read1_file.write(f"@{read_id}\n{seq1}\n+\n{qual1}\n")
        self.read2_file.write(f"@{read_id}\n{seq2}\n+\n{qual2}\n")

    def close(self) -> None:
        """Closes both FASTQ output files."""
        self.read1_file.close()
        self.read2_file.close()

    def __enter__(self):
        """Support for context manager."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close files when exiting context."""
        self.close()


import gzip
from typing import Tuple, Iterator


class FastqReader:
    """
    Iterator for reading paired-end reads from gzipped FASTQ files.

    Yields:
        Tuple[str, str, str, str, str]: (read_id, seq1, qual1, seq2, qual2)
    """

    def __init__(self, read1_filename: str, read2_filename: str, *, trim_tail: bool = False) -> None:
        """
        Initializes the FastqReader.

        Args:
            read1_filename (str): Path to R1 gzipped FASTQ file.
            read2_filename (str): Path to R2 gzipped FASTQ file.
            trim_tail : bool, optional
            If True, truncate everything in the FASTQ header line
            after the first whitespace character.  Default is False.
        """
        self.trim_tail = trim_tail
        self.read1_file = gzip.open(read1_filename, 'rt')
        self.read2_file = gzip.open(read2_filename, 'rt')

    def __iter__(self) -> Iterator[Tuple[str, str, str, str, str]]:
        return self

    def __next__(self) -> Tuple[str, str, str, str, str]:
        """
        Reads the next pair of reads.

        Returns:
            Tuple[str, str, str, str, str]: read_id, seq1, qual1, seq2, qual2

        Raises:
            StopIteration: If end of file is reached.
            RuntimeError: On unexpected error.
        """
        try:
            r1_id = next(self.read1_file).strip()
            r1_seq = next(self.read1_file).strip()
            next(self.read1_file)  # '+'
            r1_qual = next(self.read1_file).strip()

            r2_id = next(self.read2_file).strip()
            r2_seq = next(self.read2_file).strip()
            next(self.read2_file)  # '+'
            r2_qual = next(self.read2_file).strip()
            read_id = r1_id.lstrip('@')
            if self.trim_tail:
                read_id = read_id.split(" ", 1)[0]
            return read_id, r1_seq, r1_qual, r2_seq, r2_qual
        except StopIteration:
            self.close()
            raise
        except Exception as e:
            self.close()
            raise RuntimeError("Error reading FASTQ pairs") from e

    def close(self) -> None:
        """Closes both input FASTQ files."""
        self.read1_file.close()
        self.read2_file.close()

    def __enter__(self):
        """Support for context manager."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Close files when exiting context."""
        self.close()
