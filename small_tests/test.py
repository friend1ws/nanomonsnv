#! /usr/bin/env python

import threading, time

def worker(semaphore, index):
    with semaphore:
        print("worker no." + str(index) + " start!")
        time.sleep(10)
        print("worker no." + str(index) + " end!")

if __name__ == '__main__':

    semaphore = threading.Semaphore(5)
    ts = []
    for i in range(30):
        t = threading.Thread(target = worker, args = (semaphore, i,))
        t.start()
        ts.append(t)

    [t.join() for t in ts]
    print("All the threads ended!!")

