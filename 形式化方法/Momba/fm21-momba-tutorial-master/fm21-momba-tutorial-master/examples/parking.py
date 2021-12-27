# -*- coding:utf-8 -*-

from momba import jani, model
from momba.moml import expr, prop
from momba.tools import modest


FAIL_PROBABILITY = 0.1


ctx = model.Context(model.ModelType.MDP)

ctx.global_scope.declare_variable("goal", model.types.BOOL, initial_value=False)

press = ctx.create_action_type("press").create_pattern()
enter = ctx.create_action_type("enter").create_pattern()
leave = ctx.create_action_type("leave").create_pattern()


# Create the Parking Lot Automaton.

parking_lot_aut = ctx.create_automaton(name="parking-lot")

parking_lot_aut.scope.declare_variable("counter", model.types.INT, initial_value=0)

closed_loc = parking_lot_aut.create_location("closed", initial=True)
open_loc = parking_lot_aut.create_location("open")

parking_lot_aut.create_edge(
    closed_loc,
    destinations=[
        model.create_destination(
            closed_loc, probability=expr("$p", p=FAIL_PROBABILITY)
        ),
        model.create_destination(
            open_loc,
            probability=expr("1 - $p", p=FAIL_PROBABILITY),
            assignments={"counter": expr("counter + 1")},
        ),
    ],
    guard=expr("counter < 100"),
    action_pattern=press,
)

parking_lot_aut.create_edge(
    open_loc,
    destinations=[model.create_destination(open_loc)],
    guard=expr("counter >= 100"),
    action_pattern=press,
)

parking_lot_aut.create_edge(
    closed_loc,
    destinations=[
        model.create_destination(
            closed_loc, assignments={"counter": expr("counter - 1")}
        )
    ],
    action_pattern=leave,
)

parking_lot_aut.create_edge(
    open_loc,
    destinations=[model.create_destination(closed_loc)],
    action_pattern=enter,
)

parking_lot_aut.create_edge(
    open_loc,
    destinations=[model.create_destination(open_loc)],
    action_pattern=press,
)

parking_lot_aut.create_edge(
    open_loc,
    destinations=[
        model.create_destination(open_loc, assignments={"counter": expr("counter - 1")})
    ],
    action_pattern=leave,
)


# Create the Driver Automaton.

driver_aut = ctx.create_automaton(name="driver")

driver_aut.scope.declare_variable("tries", model.types.INT, initial_value=0)

press_loc = driver_aut.create_location("press", initial=True)
parked_loc = driver_aut.create_location("parked")

driver_aut.create_edge(
    press_loc,
    destinations=[
        model.create_destination(parked_loc, assignments={"goal": expr("true")})
    ],
    action_pattern=enter,
)

driver_aut.create_edge(
    press_loc,
    destinations=[
        model.create_destination(press_loc, assignments={"tries": expr("tries + 1")})
    ],
    action_pattern=press,
    guard=expr("tries < 3"),
)


# Create the automaton network.

parking_lot = parking_lot_aut.create_instance()
driver = driver_aut.create_instance()

network = ctx.create_network()

network.create_link({parking_lot: press, driver: press}, result=press)
network.create_link({parking_lot: enter, driver: enter}, result=enter)


print(jani.dump_model(network, indent=2))


checker = modest.get_checker(accept_license=True)

print(
    "Goal Probability:",
    float(
        checker.check(
            network, properties={"goal": prop("min({ Pmax(F(goal)) | initial })")}
        )["goal"]
    ),
)
