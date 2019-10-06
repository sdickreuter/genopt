import sequtils
import algorithm
import math

proc zeros*(n: int): seq[float] =
    result = newSeq[float](n)
    for i in 0..<len(result):
        result[i] = 0.0


proc ones*(n: int): seq[float] =
    result = newSeq[float](n)
    for i in 0..<len(result):
        result[i] = 1.0


proc linspace*(start, stop:float, n: int): seq[float] =
    result = newSeq[float](n)
    for i in 0..<len(result):
        result[i] = ( (float(i)/float(n) ) * (stop-start) )+start


proc arange*(start, stop, step: int): seq[int] =
    var
        n: int = int(floor( float(stop-start)/float(step)))
    result = newSeq[int](n)
    for i in 0..<len(result):
        result[i] = start+i*step


proc argsort*[T](a : T) : seq[int] {.inline.} =
  result = toSeq(0..a.len - 1)
  sort(result, proc (i, j: int): int = cmp(a[i], a[j]))


proc `[]`*[T](s: seq[T], slice: seq[int]): seq[T] {.inline.} =
    result = newSeq[T](len(slice))
    for i in 0..<result.len:
      result[i] = s[slice[i]]


proc cmp*(a, b: float): int =
    if classify(a) == fcNan:
        result = 1
    elif classify(b) == fcNan:
        result = -1
    else:
        if a > b:
            result = 1
        elif a < b:
            result = -1
        elif a == b:
            result = 0