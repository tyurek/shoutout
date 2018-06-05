from PolyCommitPed import *
import math
from helperfunctions import *
from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair

class ShoutoutRecipient:
    def __init__(self, pk, r, serverids, group, matrixdim, send_function, recv_function):
        n = len(serverids)
        t = int(math.floor((n - 1)/3))
        pc = PolyCommitPed(t=t, pk=pk, group=group, symflag=False)
        self.finished = False
        while not self.finished:
            sender, msg = recv_function()
            self.receive_msg(sender,msg)
            self.finished = True

    def receive_msg(self, sender, msg):
        if msg[0] == "payload":
            pass
