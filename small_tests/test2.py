#! /usr/bin/env python

import concurrent.futures, time, sys

def worker(index):

    print("worker no." + str(index) + " start!")
    time.sleep(4)
    print(index1)
    # print("worker no." + str(index) + " end!")
    return("worker no." + str(index) + " OK!")




if __name__ == '__main__':

    # executor = concurrent.futures.ProcessPoolExecutor(max_workers = 4)
    with concurrent.futures.ProcessPoolExecutor(max_workers = 4) as executor:

        futures = [executor.submit(worker, i) for i in range(30)]

        for future in concurrent.futures.as_completed(futures):
            try:
                ret = future.result()
            except Exception as exc:
                print(exc)
                for f in futures:
                    fret = f.cancel()
                    print(fret)
                executor.shutdown()
                sys.exit(1)
            else:
                print("OK")

            """
            if future.exception() is not None:
                print(future)
                print(future.exception())
                executor.shutdown(wait = False)
                sys.exit(1)                
            else:
                print(future.result())
            # print(future.exception())
            """

    print("All the threads ended!!")

    
