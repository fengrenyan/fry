# -*- coding:utf-8 -*-

from asciimatics.screen import Screen

# import the state space exploration engine of Momba
from momba import engine

# import the model we just defined
from . import model


def game(track: model.Track, screen: Screen, *, fail_probability: float = 0.4):
    """An interactive simulation of the jump'n'run game."""

    # TODO: Construct an automaton network by using the function
    # `model.construct_model` you defined earlier. This function
    # takes the `track` and `fail_probability` as arguments.
    network = ...

    # TODO: Construct a state space explorer from the constructed
    # automaton network. The network is a discrete time model, hence,
    # you have to use the static function `new_discrete_time` of the
    # `engine.Explorer` class. For more details see:
    # - https://momba.dev/reference/engine/#momba.engine.Explorer
    explorer = ...

    # TODO: The model has exactly one initial state. The explorer
    # allows you to query the initial states of the model. Assign
    # them to `(state,)` to extract the initial state.
    (state,) = ...

    while True:
        # TODO: We now like to print the current state to the console
        # such that the player can see the state. The code for printing
        # is provided (see below). All you have to do, is to access the
        # global variables of the state to extract the position of the
        # player. Recapitulate that the position is stored in the global
        # variables `pos_x` and `pos_y`. Hint: Have a look at the
        # properties and methods of the `State` class:
        # - https://momba.dev/reference/engine/#momba.engine.State
        #
        # Note that values have a specific type. To convert them to
        # integers use the `as_int` property:
        # - https://momba.dev/reference/engine/#momba.engine.Value.as_int
        x = ...
        y = ...

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

        # We now have to check whether there are any transitions from the
        # current state. If there are transitions, we will ask the user
        # to chose a transition, i.e., a move. Otherwise, we show a screen
        # indicating whether the player won or crashed.

        # TODO: Construct a condition for this `if` that is true if and
        # only if there are no more transitions. Hint: You again need
        # to use a property of the `State` class.
        if ...:
            # Clear the screen
            screen.clear()

            # TODO: Check whether the player won by querying the global
            # environment as you did before when extracting the position
            # of the player. This time, use `as_bool` to convert the
            # extracted value to a boolean.
            if ...:
                screen.centre("Congratulations, you won! 🏆", screen.height // 2)
            else:
                assert ...
                screen.centre("You crashed! 💥", screen.height // 2)

            screen.refresh()
            # Allow the user to quit the program by pressing any key.
            screen.wait_for_input(300)
            return

        # Wait for some user input.
        screen.wait_for_input(300)
        event = screen.get_event()
        assert event is not None

        # Allow the user to press Q to quit the game
        if event.key_code in {ord("Q"), ord("q")}:
            return

        # We now construct a dictionary mapping action labels
        # to transitions. This will allow us to later choose
        # a transition based on it's label.
        transitions = {}
        for transition in state.transitions:
            if transition.action is not None:
                # TODO: Extract the action label from the transition. You may
                # find the following pages of the documentaion helpful:
                # - https://momba.dev/reference/engine/#momba.engine.Transition
                # - https://momba.dev/reference/engine/#momba.engine.Action.action_type
                # - https://momba.dev/reference/model/actions/#momba.model.ActionType
                label = ...
                transitions[label] = transition

        # Select a transition based on the input of the user and its label.
        if event.key_code == Screen.KEY_UP:
            transition = transitions["up"]
        elif event.key_code == Screen.KEY_DOWN:
            transition = transitions["down"]
        else:
            transition = transitions["stay"]

        # TODO: Pick a successor state of the transition at random. Hint: You
        # can call `pick` on the destinations of the transition and then
        # obtain the successor state from the resulting destination:
        # - https://momba.dev/reference/engine/#momba.engine.Destination
        state = ...
