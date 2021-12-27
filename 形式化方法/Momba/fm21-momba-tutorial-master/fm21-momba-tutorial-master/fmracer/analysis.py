# -*- coding:utf-8 -*-

import typing as t

import click

from momba.moml import prop

# Import the Modest Toolset which we are using for the analysis.
from momba.tools import modest

from . import model

# Here we define the properties we would like to analyze using a model
# checker. In this case, we are interested in the maximal probability
# (`Pmax`) of eventually (`F`) winning. The function `prop` provided
# by Momba, provides a concise way to specify properties.
PROPERTIES = {"goal": prop("min({ Pmax(F(has_won)) | initial })")}


def analyze(tracks: t.Sequence[t.Tuple[str, model.Track]], resolution: int = 10):
    # Import Matplotlib, a plotting library for Python.
    import matplotlib.pyplot as plt

    # We first have to construct an instance of the model checker. This
    # will automatically download the Modest Toolset for your platform
    # and store it such that it can be invoked by Momba.
    checker = modest.get_checker(accept_license=True)

    # We now need to carry out multiple analysis on each track and plot
    # the results into a graph.
    for name, track in tracks:
        # The x axis of the plot shows the failure probability and the
        # y axis of the plot the winning probability
        x, y = [], []
        print(f"Analyzing track {name!r}...")
        # We are using `click` here to display a nice progress bar.
        with click.progressbar(range(resolution), length=resolution) as bar:
            for step in bar:
                fail_probability = step / (resolution - 1)

                # TODO: Construct an automaton network using again the
                # function `model.construct_model` you wrote earlier.
                network = ...

                # TODO: Invoke the model checker using the method `check`.
                # This method takes an automaton network as well as the
                # properties defined above (`PROPERTIES`). See:
                # - https://momba.dev/reference/analysis/#momba.analysis.Checker
                result = ...

                x.append(fail_probability)
                y.append(float(result["goal"]))
        plt.plot(x, y)

    plt.legend([name for name, _ in tracks])

    plt.ylabel("Win Probability")
    plt.xlabel("Fail Probability")

    plt.savefig("results.pdf")
