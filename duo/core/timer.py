import timeit
import duo.core.duo_exception as de

#------------------------------------------------------------
#------------------------------------------------------------
class Timer:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.start = 0.0
        self.stop = 0.0
        self.timeElapsed = 0.0

    #------------------------------------------------------------
    #------------------------------------------------------------
    def StartOrResume(self):
        self.start = timeit.default_timer()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Stop(self):
        self.stop = timeit.default_timer()
        self.timeElapsed += self.stop - self.start

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetElapsedTimeInSecond(self):
        return self.timeElapsed


#------------------------------------------------------------
#------------------------------------------------------------
class TimerManager:
    #------------------------------------------------------------
    #------------------------------------------------------------
    def __init__(self):
        self.timerList = {}

    #------------------------------------------------------------
    #------------------------------------------------------------
    def StartOrResume(self, name):
        if name not in self.timerList.keys():
            self.timerList[name] = Timer()

        self.timerList[name].StartOrResume()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Stop(self, name):
        if name not in self.timerList.keys():
            raise de.DuoException("--> Specified timer not found.")

        self.timerList[name].Stop()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def GetElapsedTimeInSecond(self, name):
        if name not in self.timerList.keys():
            raise de.DuoException("--> Specified timer not found.")

        return self.timerList[name].GetElapsedTimeInSecond()

    #------------------------------------------------------------
    #------------------------------------------------------------
    def Show(self):
        print("--> Timing result")
        for key, value in self.timerList.items():
            message = "    {0:s} : {1:f} [sec]".format(key, value.GetElapsedTimeInSecond())
            print(message)




