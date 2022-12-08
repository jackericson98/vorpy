from System.system import System, Network

file = input("Please enter the location of the system file")

mySys = System(file)

enter_network_files = input("Would you like to add a network file? (Y/N)")

if enter_network_files.lower() == "y":
    my_net_file = input("Please enter the location of the network file")
    mySys.load_net(my_net_file)
else:
    mySys.net = Network(sys=mySys, atoms=mySys.atoms)

build_net = input("\nWould you like to change ony of the build network settings before continuing?\n  Settings:\n "
                  " Surface Resolution = {}\n  Maximum Vertex Radius = {}\n  Box Size = {}\n  Solute Vertices = {}"
                  .format(mySys.net.min_dist, mySys.net.max_vert, mySys.net.box_size, mySys.net.sol_verts))


if build_net.lower() == 'y':
    changes = input("\nPlease enter the new settings as a list in order with None as placeholders "
                    "(e.g. \"[None, 5, None, None]\")")
else:
    changes = [None, None, None, None]

for i in range(len(changes)):
    if changes[i] is None or changes[i].lower() == 'none':
        changes[i] = None
    elif changes[i].lower() == 'true' or changes[i].lower() == 'false':
        changes[i] = bool(changes[i])
    else:
        changes[i] = float(changes[i])

mySys.build_network(surf_res=changes[0], max_vert=changes[1], box_size=changes[2],
                    sol_verts=changes[3])

keep_running = input("Would you like to analyze the network? (Y/N)")