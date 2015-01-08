# Awk program called in utest.sh
{
    # c = current cutsize
    # l = target cutsize
    c = $7
    s = "program=" p ": cutsize=" c " target=" t
    if (c == t) {
        print "Passed:", s;
    } else {
        print "Failed:", s;
    }
}

# EOF
