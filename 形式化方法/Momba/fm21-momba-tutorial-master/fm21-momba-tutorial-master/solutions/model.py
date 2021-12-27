# -*- coding:utf-8 -*-

from __future__ import annotations

import dataclasses as d
import typing as t

import pathlib


# import the modelling API from Momba
from momba import model
from momba.moml import expr


@d.dataclass(frozen=True)
class Cell:
    """Represents a cell on a track."""

    x: int
    y: int


@d.dataclass(frozen=True)
class Track:
    """Represents a track with a certain width, height, and obstacles."""

    width: int
    height: int
    obstacles: t.Set[Cell]

    @classmethod
    def from_ascii(cls, ascii: t.Tuple[str, ...]) -> Track:
        """Constructs a track from an ASCII art representation."""
        width = len(ascii[0])
        height = len(ascii)
        obstacles = set()
        for y, line in enumerate(ascii):
            for x, symbol in enumerate(line):
                if symbol == "x":
                    obstacles.add(Cell(x, y))
        return cls(width, height, obstacles)

    @classmethod
    def from_path(cls, path: pathlib.Path) -> Track:
        """Loads a track from a file."""
        return cls.from_ascii(tuple(path.read_text("utf-8").splitlines(keepends=False)))


def construct_model(track: Track, *, fail_probability: float = 0.4) -> model.Network:
    """
    Constructs a model for the given track.

    Here *fail_probability* is the probability that an action has no effect.
    """

    ctx = model.Context(model.ModelType.MDP)

    ctx.global_scope.declare_variable("pos_x", typ=model.types.INT, initial_value=0)
    ctx.global_scope.declare_variable("pos_y", typ=model.types.INT, initial_value=0)

    ctx.global_scope.declare_variable(
        "has_won", typ=model.types.BOOL, initial_value=False
    )
    ctx.global_scope.declare_variable(
        "has_crashed", typ=model.types.BOOL, initial_value=False
    )

    up_action = ctx.create_action_type("up")
    down_action = ctx.create_action_type("down")
    stay_action = ctx.create_action_type("stay")

    automaton = ctx.create_automaton(name="Environment")

    location = automaton.create_location("ready", initial=True)

    can_move = expr("not has_won and not has_crashed")

    def has_won(x: model.Expression) -> model.Expression:
        """Constructs an expression for whether the player won."""
        return expr("$x >= $width", x=x, width=track.width)

    def has_crashed(x: model.Expression, y: model.Expression) -> model.Expression:
        """Constructs an expression for whether the player crashed."""
        out_of_bounds = expr("$y >= $height or $y < 0", y=y, height=track.height)
        on_obstacle = model.expressions.logic_any(
            *(
                expr(
                    "$x == $obstacle_x and $y == $obstacle_y",
                    x=x,
                    y=y,
                    obstacle_x=obstacle.x,
                    obstacle_y=obstacle.y,
                )
                for obstacle in track.obstacles
            )
        )
        return model.expressions.logic_or(out_of_bounds, on_obstacle)

    moves = {up_action: -1, down_action: 1, stay_action: 0}

    for action_type, delta in moves.items():
        new_pos_x = expr("pos_x + 1")
        new_pos_y = expr("pos_y + $delta", delta=delta)

        success_probability = expr(
            "1 - $fail_probability", fail_probability=fail_probability
        )

        success_destination = model.create_destination(
            location,
            probability=success_probability,
            assignments={
                "pos_x": new_pos_x,
                "pos_y": new_pos_y,
                "has_won": has_won(new_pos_x),
                "has_crashed": has_crashed(new_pos_x, new_pos_y),
            },
        )

        fail_destination = model.create_destination(
            location,
            probability=expr("$p", p=fail_probability),
            assignments={
                "pos_x": new_pos_x,
                "has_won": has_won(new_pos_x),
                "has_crashed": has_crashed(new_pos_x, expr("pos_y")),
            },
        )

        automaton.create_edge(
            location,
            action_pattern=action_type.create_pattern(),
            guard=can_move,
            destinations=[success_destination, fail_destination],
        )

    network = ctx.create_network()

    instance = automaton.create_instance()
    network.add_instance(instance)

    for action_type in moves.keys():
        network.create_link(
            {instance: action_type.create_pattern()},
            result=action_type.create_pattern(),
        )

    return network
