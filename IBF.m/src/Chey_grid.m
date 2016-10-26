function grid = Chey_grid(NG)

grid = (cos((NG-1:-1:0)/(NG-1)*pi)+1)/2;
grid = grid(:);

end