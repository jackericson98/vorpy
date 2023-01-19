from System.system import *
from Visualize.commands.build import *
from Visualize.commands.export import *
from Visualize.commands.load import *
from Visualize.commands.set import *


def check_input():
    """
    Main function that is looped. Checks the inputs and runs the correct functions
    :return:
    """
    global sys
    # Set up the prompt
    usr_npt = input("vorpy >>>   ")
    # Split it up by the spaces
    usr_npt = usr_npt.split()
    # Check to see if the initial input is a command
    print(usr_npt[0].lower())
    if len(usr_npt) == 0 or usr_npt[0].lower() not in my_commands:
        invalid_input(usr_npt)
        return True

    ########################## Commands  ################################################

    # Check if the user's input is in loads
    if usr_npt[0].lower() in load_cmds:
        print(usr_npt)
        my_sys = load(sys=sys, usr_npt=usr_npt)
        if my_sys is not None:
            sys = my_sys
            return sys
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in set_cmds:
        sett(sys=sys, usr_npt=usr_npt)
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in build_cmds:
        build(sys=sys)
    # Check if the user's input is in exports
    elif usr_npt[0].lower() in export_cmds:
        export(sys=sys, usr_npt=usr_npt)
    # Check if the user's input is in shows
    elif usr_npt[0].lower() in show_cmds:
        show(sys=sys, usr_npt=usr_npt)
    # Check if the user wants help
    elif usr_npt[0].lower() in helps:
        help_()
    # Check to see if the user's input includes quits
    elif usr_npt[0].lower() in quits:
        if are_you_sure():
            return False

    # Unless the user quits, keep running the program
    return True


if __name__ == '__main__':
    # Welcome introduction
    print("Welcome to vorpy. For assistance type \'h\'. To quit type \'q\'")
    # Create the system
    sys = System()
    # Run the program
    while True:
        # create_header(mySys)
        running = check_input()
        if type(running) is not bool:
            sys = running
        else:
            break
