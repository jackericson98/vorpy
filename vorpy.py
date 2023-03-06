import sys

from Visualize.commands.build import *
from Visualize.commands.export import *
from Visualize.commands.load import *
from Visualize.commands.set import *
from Visualize.commands.argv_commands import interpret_argvs


def check_input():
    """
    Main function that is looped. Checks the inputs and runs the correct functions
    :return:
    """
    global my_sys
    # Set up the prompt
    usr_npt = input("vorpy >>>   ")
    # Split it up by the spaces
    usr_npt = usr_npt.split()
    # Check to see if the initial input is a command
    if len(usr_npt) == 0 or usr_npt[0] not in my_commands:
        invalid_input(usr_npt)
        return True

    ########################## Commands  ################################################

    # Check if the user's input is in loads
    if usr_npt[0].lower() in load_cmds:
        my_sys = load(sys=my_sys, usr_npt=usr_npt)
        if my_sys is not None:
            sys = my_sys
            return sys
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in set_cmds:
        sett(sys=my_sys, usr_npt=usr_npt)
    # Check if the user's input is in builds
    elif usr_npt[0].lower() in build_cmds:
        if len(usr_npt) == 1:
            usr_npt = None
        build(sys=my_sys, usr_npt=usr_npt)
    # Check if the user's input is in exports
    elif usr_npt[0].lower() in export_cmds:
        export(sys=my_sys, usr_npt=usr_npt)
    # Check if the user's input is in shows
    elif usr_npt[0].lower() in show_cmds:
        show(sys=my_sys, usr_npt=usr_npt)
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
    # Create the system
    my_sys = System(root_dir=os.getcwd())
    # Check for arguments from the user
    if len(sys.argv) > 1:
        interpret_argvs(my_sys)
    else:
        # Welcome introduction
        print("Welcome to vorpy. For assistance type \'h\'. To quit type \'q\'")
        # Set up the running variable
        running = True
        # Run the program
        while running:
            # create_header(mySys)
            running = check_input()
            if type(running) is not bool:
                my_sys = running
                running = True
