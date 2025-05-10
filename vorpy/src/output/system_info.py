

def export_sys_info(sys):
    """
    Exports system information to a text file.
    
    This function creates a detailed information file about the system, including:
    - System name and network information
    - Chain information including atom counts, residue counts, volumes, and surface areas
    - Group information including atom counts, residue counts, chain counts, volumes, and surface areas
    
    The information is written to a file named {system_name}_info.txt in the system's output directory.
    
    :param sys: System object containing the information to be exported
    :return: None
    """
    # Open the file
    with open(sys.name + "_info.txt", 'w') as info:
        # Write the header
        info.write(sys.name + " Network")
        # Write the chain header
        info.write("\n\n++++++++++++++++++++++++  Chains  +++++++++++++++++++++++++++++++\n\n")
        # Go through the chains in the system
        if sys.chains is not None:
            for chain in sys.chains:
                # Write the chain header
                info.write("Chain {} - {} atoms, {} residues\n\n".format(chain.name, len(chain.atoms), len(chain.residues)))
                # Quick check to see if the chain has been calculated
                if chain.vol is not None and chain.vol < 0:
                    # Write the chain information
                    info.write("  Volume = {}, Surface Area = {}\n\n\n".format(chain.vol, chain.sa))
        # Draw a separating line
        info.write("\n\n++++++++++++++++++++++++  Groups  +++++++++++++++++++++++++++++++\n\n")
        for group in sys.groups:
            # Write the group header
            info.write("Group {} - {} atoms, {} residues, {} chains\n\n".format(group.name, len(group.atms), len(group.rsds), len(group.chns)))
            # Write the group info
            info.write("  Volume = {}, Surface Area = {}\n\n\n".format(group.vol, group.sa))