import sqlite3
import tempfile
from pathlib import Path
from types import TracebackType
from typing import Iterator, Literal


CarrierKeyOrder = Literal["feature", "sample"]


class MinimumDistanceStore:
    """Disk-backed minimum reduction for carrier keys."""

    def __init__(self, directory: Path):
        temporary = tempfile.NamedTemporaryFile(
            prefix="carrier-minima-",
            suffix=".sqlite3",
            dir=directory,
            delete=False,
        )
        temporary.close()
        self.path = Path(temporary.name)
        self.connection = sqlite3.connect(self.path)
        self.connection.execute("PRAGMA journal_mode=OFF")
        self.connection.execute("PRAGMA synchronous=OFF")
        self.connection.execute("PRAGMA temp_store=FILE")
        self.connection.execute("PRAGMA cache_size=-4096")
        self.connection.execute(
            """
            CREATE TABLE carrier_minima (
                feature_id TEXT NOT NULL,
                sample_id TEXT NOT NULL,
                ac_class TEXT NOT NULL,
                minimum_distance_bp INTEGER NOT NULL CHECK (minimum_distance_bp >= 0),
                PRIMARY KEY (feature_id, sample_id, ac_class)
            ) WITHOUT ROWID
            """
        )

    def __enter__(self) -> "MinimumDistanceStore":
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def close(self) -> None:
        self.connection.close()
        self.path.unlink(missing_ok=True)

    def upsert(
        self,
        sample_id: str,
        feature_id: str,
        ac_class: str,
        minimum_distance_bp: int,
    ) -> None:
        self.connection.execute(
            """
            INSERT INTO carrier_minima (
                feature_id, sample_id, ac_class, minimum_distance_bp
            ) VALUES (?, ?, ?, ?)
            ON CONFLICT (feature_id, sample_id, ac_class) DO UPDATE SET
                minimum_distance_bp = min(
                    carrier_minima.minimum_distance_bp,
                    excluded.minimum_distance_bp
                )
            """,
            (feature_id, sample_id, ac_class, minimum_distance_bp),
        )

    def count(self) -> int:
        row = self.connection.execute("SELECT count(*) FROM carrier_minima").fetchone()
        assert row is not None
        return int(row[0])

    def distinct_feature_count(self) -> int:
        row = self.connection.execute(
            "SELECT count(DISTINCT feature_id) FROM carrier_minima"
        ).fetchone()
        assert row is not None
        return int(row[0])

    def iter_feature(self, feature_id: str) -> Iterator[tuple[str, str, int]]:
        cursor = self.connection.execute(
            """
            SELECT sample_id, ac_class, minimum_distance_bp
            FROM carrier_minima
            WHERE feature_id = ?
            ORDER BY sample_id, ac_class
            """,
            (feature_id,),
        )
        for sample_id, ac_class, distance in cursor:
            yield str(sample_id), str(ac_class), int(distance)

    def write_tsv(self, path: Path, key_order: CarrierKeyOrder) -> None:
        if key_order == "feature":
            order_by = "feature_id, sample_id, ac_class"
        elif key_order == "sample":
            order_by = "sample_id, feature_id, ac_class"
        else:
            raise ValueError(f"Unsupported carrier key order: {key_order}")

        with path.open("w", encoding="utf-8") as handle:
            handle.write("sample_id\tfeature_id\tac_class\tminimum_distance_bp\n")
            cursor = self.connection.execute(
                """
                SELECT sample_id, feature_id, ac_class, minimum_distance_bp
                FROM carrier_minima
                ORDER BY """
                + order_by
            )
            for sample_id, feature_id, ac_class, distance in cursor:
                handle.write(f"{sample_id}\t{feature_id}\t{ac_class}\t{distance}\n")
