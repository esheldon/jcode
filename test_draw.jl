using gmix
using Winston
using Colors

#function test_render(model; T=16.0, g1=0.2, g2=0.3, show=false)


    model=gmix.GAUSS
    T=16.0
    g1=0.2
    g2=0.3
    show=false

    flux=100.0
    sigma=sqrt(T/2.0)
    dim=ceil(Int, 2.0*5.0*sigma)

    dims = (dim,dim)

    pars = [dim/2.0, dim/2.0, g1, g2, T, flux]

    gm = gmix.make_simple(model, pars)
    println("gm:")
    println(gm)

    im = gmix.make_image(gm, dims)

    if show
        cmap = reverse(Colors.colormap("grays"))
        Winston.colormap(cmap)
        implt=Winston.imagesc(im)
        setattr(implt, aspect_ratio=1.0)

        Winston.display(implt)
        junk = readline(STDIN)
    end

#end
