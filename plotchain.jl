module plotchain

using Winston

type Grid
    nplot::Int
    nrow::Int
    ncol::Int
    nplot_tot::Int

    function Grid(nplot)
        self=new(0,0,0)

        self.nplot = floor(Int, nplot)
        

        # first check some special cases
        if self.nplot==8
            self.nrow, self.ncol = 2,4
        else

            sq=floor(Int, sqrt(self.nplot))

            if self.nplot==sq*sq
                self.nrow, self.ncol = sq,sq
            elseif nplot <= sq*(sq+1)
                self.nrow, self.ncol = sq,sq+1
            else
                self.nrow, self.ncol = sq+1,sq+1
            end
        end
        self.nplot_tot=self.nrow*self.ncol

        self
    end
end

function get_rowcol(self::Grid, index)
    index = floor(Int, index)

    if index > self.nplot_tot
        throw(error("index too large $(index) > $(self.nplot_tot)"))
    end

    row0 = div(index-1,self.ncol)
    col0 =  (index-1) % self.ncol

    return 1+row0,1+col0
end

function plot(chain; nbin=50, show=false)

    ntrials, ndim = size(chain)

    grid = Grid(ndim)
    tab = Table(grid.nrow, grid.ncol)

    for i=1:ndim
        p = FramedPlot(xlabel="p$(i)",
                       ylabel="N")

        add(p, Histogram(hist(chain[:,i], nbin)...))

        row,col = get_rowcol(grid, i)
        tab[row,col] = p
       
    end

    if show
        display(tab)
        print("hit a key: ")
        junk = readline(STDIN)
    end
    tab
end

end # module
