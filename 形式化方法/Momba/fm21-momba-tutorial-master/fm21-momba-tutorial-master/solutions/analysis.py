# -*- coding:utf-8 -*-

import typing as t

import click


from momba.moml import prop
from momba.tools import modest

from . import model


PROPERTIES = {"goal": prop("min({ Pmax(F(has_won)) | initial })")}


def analyze(tracks: t.Sequence[t.Tuple[str, model.Track]], resolution: int = 10):
    import matplotlib.pyplot as plt

    checker = modest.get_checker(accept_license=True)

    for name, track in tracks:
        x = []
        y = []
        print(f"Analyzing track {name!r}...")
        with click.progressbar(range(resolution), length=resolution) as bar:
            for step in bar:
                fail_probability = step / (resolution - 1)

                network = model.construct_model(
                    track, fail_probability=fail_probability
                )
                result = checker.check(network, properties=PROPERTIES)

                x.append(fail_probability)
                y.append(float(result["goal"]))
        plt.plot(x, y)

    plt.legend([name for name, _ in tracks])

    plt.ylabel("Win Probability")
    plt.xlabel("Fail Probability")

    plt.savefig("results.pdf")
