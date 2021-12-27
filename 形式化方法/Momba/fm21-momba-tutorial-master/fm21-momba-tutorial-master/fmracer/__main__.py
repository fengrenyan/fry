# -*- coding:utf-8 -*-

"""Command line tool."""

import typing as t

import functools
import pathlib
import warnings

import click

from asciimatics.screen import Screen

from . import analysis, model, interactive

warnings.simplefilter("ignore")


@click.group()
def main():
    pass


@main.command()
@click.argument("tracks", metavar="TRACK", type=pathlib.Path, required=True, nargs=-1)
@click.option("--resolution", type=int, default=10)
def analyze(tracks: t.Sequence[pathlib.Path], resolution: int):
    """Analyze the winning probabilities for the given tracks."""
    analysis.analyze(
        [(track.stem, model.Track.from_path(track)) for track in tracks], resolution
    )


@main.command()
@click.argument("track", type=pathlib.Path)
@click.option("--fail-probability", type=float, default=0.4)
def race(track: pathlib.Path, fail_probability: float):
    """Start the interactive simulation of the game."""
    Screen.wrapper(
        functools.partial(
            interactive.game,
            model.Track.from_path(track),
            fail_probability=fail_probability,
        )
    )


if __name__ == "__main__":
    main()
