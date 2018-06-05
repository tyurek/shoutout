from ShoutoutDealer import *
from ShoutoutRecipient import *
import gevent
from gevent import Greenlet
from gevent.queue import Queue
import random

def simple_router(participantids, maxdelay=0.01, seed=None):
    """Builds a set of connected channels, with random delay
    @return (receives, sends)
    """
    rnd = random.Random(seed)
    #if seed is not None: print 'ROUTER SEED: %f' % (seed,)
    
    queues = {}
    for _ in participantids:
        queues[_] = Queue()
    
    def makeSend(i):
        def _send(j, o):
            delay = rnd.random() * maxdelay
            #print 'SEND %8s [%2d -> %2d] %.2f' % (o[0], i, j, delay)
            gevent.spawn_later(delay, queues[j].put, (i,o))
            #queues[j].put((i, o))
        return _send
    
    def makeRecv(j):
        def _recv():
            (i,o) = queues[j].get()
            #print 'RECV %8s [%2d -> %2d]' % (o[0], i, j)
            return (i,o)
        return _recv
    
    sends = {}
    receives = {}
    for i in participantids:
        sends[i] = makeSend(i)
        receives[i] = makeRecv(i)    
    return (sends, receives)

def gen_pubkey(t, group):
    alpha = group.random(ZR)
    g = group.random(G1)
    ghat = group.random(G2)
    h = group.random(G1)
    pk = []
    for i in range(t+1):
        pk.append(g**(alpha**i))
    for i in range(2):
        pk.append(ghat**(alpha**i))
    for i in range(t+1):
        pk.append(h**(alpha**i))
    return pk


def main():
    group = PairingGroup('BN254')
    serverids = [1,2,3,4,5,6,7]
    secret = [42,43]
    r = len(secret)
    n = len(serverids)
    t = int(math.floor((n - 1)/3))
    pk = gen_pubkey(t, group)

    dealerid = 42
    sends, recvs = simple_router(serverids + [dealerid])
    threads = []
    thread = Greenlet(ShoutoutDealer, pk=pk, secret=secret, position=4, serverids=serverids, group=group, matrixdim=6, send_function=sends[dealerid], recv_function=recvs[dealerid])
    thread.start()
    threads.append(thread)
    for pid in serverids:
        thread = Greenlet(ShoutoutRecipient, pk=pk, r=r, serverids=serverids, group=group, matrixdim=6, send_function=sends[pid], recv_function=recvs[pid])
        thread.start()
        threads.append(thread)
    gevent.joinall(threads)


if __name__ == "__main__":
    debug = True
    main()
