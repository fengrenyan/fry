# -*- coding:utf-8 -*-

from asciimatics.screen import Screen

# import the state space exploration engine of Momba
from momba import engine

from . import model


def game(track: model.Track, screen: Screen, *, fail_probability: float = 0.4):
    network = model.construct_model(track, fail_probability=fail_probability)

    explorer = engine.Explorer.new_discrete_time(network)

    (state,) = explorer.initial_states

    while True:
        x = state.global_env["pos_x"].as_int
        y = state.global_env["pos_y"].as_int
        player = model.Cell(x, y)

        # Prints the current state of the game to the console.
        for y in range(track.height):
            for x in range(track.width):
                cell = model.Cell(x, y)
                if cell in track.obstacles:
                    screen.print_at(
                        "█", x, y, Screen.COLOUR_RED, bg=Screen.COLOUR_WHITE
                    )
                elif cell == player:
                    screen.print_at(
                        "x", x, y, Screen.COLOUR_BLUE, bg=Screen.COLOUR_WHITE
                    )
                else:
                    screen.print_at(" ", x, y, bg=Screen.COLOUR_WHITE)
        screen.refresh()

        if not state.transitions:
            screen.clear()
            if state.global_env["has_won"].as_bool:
                screen.centre("Congratulations, you won! 🏆", screen.height // 2)
            else:
                assert state.global_env["has_crashed"].as_bool
                screen.centre("You crashed! 💥", screen.height // 2)
            screen.refresh()
            screen.wait_for_input(300)
            return

        screen.wait_for_input(300)
        event = screen.get_event()
        assert event is not None

        if event.key_code in {ord("Q"), ord("q")}:
            return

        transitions = {}
        for transition in state.transitions:
            if transition.action is not None:
                label = transition.action.action_type.label
                transitions[label] = transition

        if event.key_code == Screen.KEY_UP:
            transition = transitions["up"]
        elif event.key_code == Screen.KEY_DOWN:
            transition = transitions["down"]
        else:
            transition = transitions["stay"]

        state = transition.destinations.pick().state
