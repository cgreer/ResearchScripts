from cgAutoCast import autocast

@autocast
def addStuff(one, two, three, four):
        print one + two + three
        print type(four)
addStuff('1', '2', '3', 'True')        
