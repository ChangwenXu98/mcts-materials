"""
Drive the superhydride MCTS campaign, writing a snapshot after every chunk.

Runs the search in small chunks rather than one call so that results are on
disk continuously: the QE evaluator's CSV cache is already appended per
candidate, and this adds a search-level snapshot (summary, ranking,
convergence, tree) so progress can be read and plotted while the job is still
running.

    python run_campaign.py <config.yaml> <output_dir> [chunk_size]
"""

import asyncio
import json
import logging
import sys
import time
from pathlib import Path

from mcts_framework.cli.builders import build_mcts
from mcts_framework.cli.results import save_results
from mcts_framework.core.config import Config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("campaign")


def write_progress(output_dir: Path, mcts, started: float, chunk: int, done: bool) -> None:
    """A small JSON the monitoring side can poll without parsing the tree."""
    evaluated = len(mcts.property_evaluator.cache)
    best = mcts.get_best_materials(n=1)
    payload = {
        "elapsed_s": round(time.time() - started, 1),
        "chunks_done": chunk,
        "iterations": len(mcts.reward_history),
        "unique_materials": len(mcts.visited_materials),
        "evaluated": evaluated,
        "best_reward": mcts.best_reward,
        "best_formula": best[0].material.get_formula() if best else None,
        "root_terminated": bool(mcts.root.terminated),
        "exhausted": bool(mcts.terminated),
        "finished": done,
    }
    (output_dir / "progress.json").write_text(json.dumps(payload, indent=2))
    logger.info(
        "snapshot: %d iters, %d unique, %d evaluated, best %s at reward %.4f",
        payload["iterations"], payload["unique_materials"], evaluated,
        payload["best_formula"], payload["best_reward"],
    )


async def main() -> int:
    config_path, output_dir = sys.argv[1], Path(sys.argv[2])
    chunk_size = int(sys.argv[3]) if len(sys.argv) > 3 else 5

    config = Config.from_yaml(config_path)
    output_dir.mkdir(parents=True, exist_ok=True)

    mcts = build_mcts(config)
    logger.info(
        "root %s | palette %s | evaluator %s",
        mcts.root.material.get_identifier(),
        config.superhydride.host_palette,
        config.superhydride.evaluator,
    )

    preflight = getattr(mcts.property_evaluator, "preflight", None)
    if preflight is not None:
        found = preflight()
        logger.info("QE preflight: %s", found)
        if not all(found.values()):
            logger.error("QE binaries missing; every reward would be 0.0. Aborting.")
            return 1

    started = time.time()
    total = config.mcts.iterations
    chunk = 0
    while len(mcts.reward_history) < total:
        remaining = total - len(mcts.reward_history)
        await mcts.run(iterations=min(chunk_size, remaining))
        chunk += 1
        save_results(mcts, str(output_dir), top_n=50, config=config, save_tree=True)
        write_progress(output_dir, mcts, started, chunk, done=False)
        if mcts.terminated or (mcts.root.terminated and config.mcts.search_mode == "fast"):
            logger.info("search stopped early (exhausted=%s)", mcts.terminated)
            break

    save_results(mcts, str(output_dir), top_n=50, config=config, save_tree=True)
    write_progress(output_dir, mcts, started, chunk, done=True)
    logger.info("campaign complete in %.1f s", time.time() - started)
    return 0


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))
