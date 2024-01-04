function result = TupleMultiply(tuple1, tuple2)
    result = [max(tuple1(1) + tuple2(1), tuple1(2) + tuple2(2)), max(tuple1(1) + tuple2(2), tuple1(2) + tuple2(1))];
end