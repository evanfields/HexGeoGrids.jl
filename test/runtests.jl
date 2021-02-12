using HexGeoGrids
HGG = HexGeoGrids
Hexagons = HexGeoGrids.Hexagons
using Test
using Random

# controls for randomization and default tests systems
Random.seed!(47)
n_points = 1000
systems = HexSystem.(
    [0, -5, 45, -180, 180, 47.47], # longitude
    [0, 1, 45, 89, 3, -89], # latitude
    [5, 500, 5000, 1, .6, 50_000], # size
)

# haversine distances useful for testing
const EARTH_RAD_EQ = 6378137.0
const EARTH_RAD_POLAR = 6356752.3
function haversine(lon1, lat1, lon2, lat2)
	avglat = (lat1 + lat2) * pi / 360
	rad = sqrt(
	    ((EARTH_RAD_EQ^2 * cos(avglat))^2 + (EARTH_RAD_POLAR^2 * sin(avglat))^2) /
	    ((EARTH_RAD_EQ * cos(avglat))^2 + (EARTH_RAD_POLAR * sin(avglat))^2)
    )

    return 2 * rad * asin(sqrt(sind((lat2 - lat1) / 2)^2 +
    cosd(lat1) * cosd(lat2) * sind((lon2 - lon1) / 2)^2))
end

@testset "Constructors" begin
    # expected behavior
    @test_nowarn HexSystem(1, 2, 3)
    @test_nowarn HexSystem(-1, Int32(2), UInt16(3))
    @test_nowarn HexSystem(1.1, 2, 3) # outer constructor
    @test_nowarn HexSystem(0.5, 89.99, 1.02)
    @test_throws DomainError HexSystem(0, -90, 11) # poles not allowed
    @test_throws DomainError HexSystem(0, -90, 0) # size must be >= 1
end

@testset "Lon-lat to hexagon" begin
    lons = 360 * rand(n_points) .- 180
    lats = 180 * rand(n_points) .- 90
    for sys in systems, (lon, lat) in zip(lons, lats)
        if abs(lon - sys.center.lon) >= 40
            continue # can't encode this lat-lon in this hex system
        end
        # The center of the hexagon with this point should be near the point
        hex = HexCell(lon, lat, sys)
        center_point = center(hex)
        @test haversine(center_point..., lon, lat) <= sys.hex_size * 1.05
        # return to hex form should recover the same hexagon
        hex_of_center = HexCell(center_point..., sys)
        @test hex_of_center == hex
    end
end

@testset "Hex to hash" begin
    for sys in systems, _ in 1:n_points
        q, r = rand(-1000:1000, 2)
        cell = HexCell(sys, Hexagons.HexFlatTop(Hexagons.CoordAxial(q, r)))
        @test HexCell(index(cell)) == cell
    end
end

@testset "Utils" begin
    @test HGG._shift_needed(2) == 1
    @test HGG._shift_needed(17) == -2
    @test HGG._shift_needed(18) == 3
end