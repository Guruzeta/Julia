function summ(x,y=2)
    z=3
    x+y+z
end

using Debugger

#@enter summ(1)



# a new way to write functions
summ= function (a,b,c)
    a+b+c
end


#debugging a for loop
f3 = function()
    sum=0
    for i=1:20
        if i==2
            @bp
        end
        sum=sum+i^2
    end
end
@run f3()
