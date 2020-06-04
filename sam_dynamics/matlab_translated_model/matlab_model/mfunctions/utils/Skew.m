function S = Skew(lambda)

if (length(lambda) ~= 3),error('lambda-vector must have dimension 3 !');end

l1 = lambda(1);
l2 = lambda(2);
l3 = lambda(3);

S = [0 -l3 l2;
     l3 0 -l1;
     -l2 l1 0];

end