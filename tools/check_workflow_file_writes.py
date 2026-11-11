#!/usr/bin/env python3
"""Reject file-writing expressions that Terra would evaluate in a workflow.

Task expressions are allowed. The WDL parser supplies the scope and expression
tree, so comments, quoted text, scatters, and call inputs are handled correctly.
Requires miniwdl (the ``WDL`` Python package).
"""

import argparse
from collections.abc import Iterator
from pathlib import Path
import sys

import WDL


def workflow_file_writes(document: WDL.Tree.Document) -> Iterator[WDL.Expr.Apply]:
    """Find file-writing calls in this document and its imported workflows."""
    visited: set[int] = set()

    def expressions(node):
        if isinstance(node, WDL.Expr.Apply) and node.function_name.startswith("write_"):
            yield node
        for child in node.children:
            yield from expressions(child)

    def workflows(current):
        if id(current) in visited:
            return
        visited.add(id(current))
        if current.workflow is not None:
            # A call's children are its input expressions, not its task body.
            yield from expressions(current.workflow)
        for imported in current.imports:
            yield from workflows(imported.doc)

    yield from workflows(document)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wdl", type=Path, nargs="+", help="WDL file to check")
    arguments = parser.parse_args(argv)
    failures = 0
    for path in arguments.wdl:
        document = WDL.load(str(path))
        for expression in workflow_file_writes(document):
            position = expression.pos
            print(
                f"{position.uri}:{position.line}:{position.column}: "
                f"{expression.function_name} is not allowed at workflow scope; "
                "create the file within a task",
                file=sys.stderr,
            )
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
