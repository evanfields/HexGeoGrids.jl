# dummy line so that VS code will get the scope right
module HexGeoGrids

export
    HexCell,
    HexSystem,
    center,
    vertices,
    index

using Geodesy:
    LatLon,
    LLA,
    LLAfromUTM,
    UTM,
    UTMfromLLA,
    UTMZfromLLA,
    utm_zone,
    wgs84

using GeoInterface
import GeoInterface: Polygon
using Hexagons

struct UTMZone
    zone::Int
    isnorth::Bool
end

struct HexSystem
    center::LLA{Float64}
    utm_zone::UTMZone
    hex_size::UInt16
    
    function HexSystem(center_lon::Integer, center_lat::Integer, size::Integer)
        if !(-90 < center_lat < 90)
            throw(DomainError(
                center_lat,
                "May not construct HexSystems centered at the poles."
            ))
        end
        if size < 1
            throw(DomainError(size, "Size must be positive."))
        end
        center = LLA(center_lat, center_lon, 0.0)
        return new(
            center,
            UTMZone(utm_zone(center)...),
            UInt16(size)
        )
    end
end

function HexSystem(lon::Real, lat::Real, size::Real)
    # don't allow centering at the poles
    lon = round(Int, lon)
    lat = round(Int, clamp(lat, -89, 89))
    size = round(Int, size)
    return HexSystem(lon, lat, size)
end

struct HexCell
    system::HexSystem
    hex::HexFlatTop{CoordAxial{Int64}}
end 

###
# Conversions lonlat ⇔ hexagon
###



"""To ensure HexCells are nearly north aligned everywhere, longitude-shift lon-lat
points so that UTM-LLA conversions happen as if the HexSystem center is on a line
of longitude 3 mod 6. E.g. if a HexSystem has a center at longitude 7 degrees, shift
longitudes by 2 degrees so that the coordinate conversions happen near 9 degrees
of longitude."""
_shift_needed(lon) = 3 - mod(lon, 6)
_shift_needed(hs::HexSystem) = _shift_needed(hs.center.lon)

_shift_lon(lla::LLA, dlon) = LLA(lla.lat, lla.lon + dlon, lla.alt)


"""Map a latitude and longitude into cartesian coordinates for a given
HexSystem. In this coordinate system, the HexSystem's center is at the
origin."""
function _cart_coords(lon, lat, hs::HexSystem)
    shift = _shift_needed(hs)
    lla_point = _shift_lon(LLA(lat, lon, 0.0), shift)
    converter = UTMfromLLA(hs.utm_zone.zone, hs.utm_zone.isnorth, wgs84)
    utmz_center = converter(_shift_lon(hs.center, shift))
    utmz_point = converter(lla_point)
    return utmz_point.x - utmz_center.x, utmz_point.y - utmz_center.y
end

"""Map lon-lat coordinates (in degrees) to a HexCell."""
function HexCell(lon::Real, lat::Real, hs::HexSystem)
    coords = _cart_coords(lon, lat, hs) ./ hs.hex_size
    return HexCell(hs, hex_containing(HexFlatTop, coords...))
end

"""Map an (x,y) point in the Cartesian plane of Hexagons.jl to a lon-lat pair
by way of the UTM projection specified by a hex system. Return a tuple (lon, lat)."""
function _hex_cartesian_to_lonlat(cart_x, cart_y, hs::HexSystem)
    shift = _shift_needed(hs)
    latlon_to_utm = UTMfromLLA(hs.utm_zone.zone, hs.utm_zone.isnorth, wgs84)
    utm_center = latlon_to_utm(_shift_lon(hs.center, shift))
    utm_point = UTM(cart_x + utm_center.x, cart_y + utm_center.y)
    utm_to_latlon = LLAfromUTM(hs.utm_zone.zone, hs.utm_zone.isnorth, wgs84)
    lla_point = _shift_lon(utm_to_latlon(utm_point), -shift)
    return lla_point.lon, lla_point.lat
end

"""Return the vertices of a HexCell as a list of (lon, lat) tuples."""
function vertices(cell::HexCell)
    system = cell.system
    size = system.hex_size
    return [
        _hex_cartesian_to_lonlat((size .* vert)..., system)
        for vert in Hexagons.vertices(cell.hex)
    ]
end

"""Return the center of a HexCell as a tuple (lon, lat)."""
function center(cell::HexCell)
    center_cartesian = Hexagons.center(cell.hex) .* cell.system.hex_size
    return _hex_cartesian_to_lonlat(center_cartesian..., cell.system)
end

###
# Conversions hexagon ⇔ index
###

"""Convert integer longitude in [-180, 180] and integer latitude in [-90, 90] to a
single integer in [0, 65340]."""
function _lonlat_to_int(lon::Integer, lat::Integer)
    -180 <= lon <= 180 || throw(DomainError(lon, "Lon of $lon not in [-180, 180]."))
    -90 <= lat <= 90 || throw(DomainError(lat, "Lat of $lat not in [-90, 90]"))
    return (lat + 90) + 181 * (lon + 180)
end
"""Convert an integer in [0, 65340] to a lon-lat pair. Inverse of `_lonlat_to_int`."""
function _int_to_lonlat(x::Integer)
    # 65_340 is 361 * 181 - 1, where 361 is number of legal lons, 181 legal lats
    0 <= x <= 65340 || throw(DomainError(x, "x not in [0, 65_340]"))
    lon_shifted, lat_shifted = divrem(x, 181)
    return lon_shifted - 180, lat_shifted - 90
end

"""Compute the 8 digit (base 16) prefix corresponding to a HexSystem."""
function _system_to_prefix(hs::HexSystem)
    lonlat_int = _lonlat_to_int(round(Int, hs.center.lon), round(Int, hs.center.lat))
    return string(lonlat_int, base = 16, pad = 4) * string(hs.hex_size, base = 16, pad = 4)
end

"""Construct a HexSystem from an 8-digit (base 16) hash prefix."""
function _prefix_to_system(str)
    p = ascii(str)
    if length(p) != 8
        throw(DomainError("Length of argument to `_prefix_to_system` must be 8."))
    end
    lon, lat = parse(Int, p[1:4], base = 16) |> _int_to_lonlat
    hex_size = parse(UInt16, p[5:8], base = 16)
    return HexSystem(lon, lat, hex_size)
end

"""Create an 8-digit (base 16) hash body from a hexagon
(i.e. the `.hex` field of a `HexCell`)."""
function _hex_coords_to_hash(hex)
    h_axial = hex.coords
    # shift Int16 to UInt16 so string form won't have a leading '-'
    q_offset, r_offset = ((h_axial.q, h_axial.r) .- typemin(Int16)) .|> UInt16
    return string(q_offset, base = 16, pad = 4) * string(r_offset, base = 16, pad = 4)
end

"""Construct hexagon coordinates from an 8-digit (base 16) hash body."""
function _hash_to_hex_coords(str)
    q_offset = parse(Int, str[1:4], base = 16)
    r_offset = parse(Int, str[5:8], base = 16)
    return CoordAxial(q_offset + typemin(Int16), r_offset + typemin(Int16))
end

"""Encode a HexCell as a 17 character `prefix:body` index string."""
function index(hc::HexCell)
    prefix = _system_to_prefix(hc.system)
    body = _hex_coords_to_hash(hc.hex)
    return prefix * ":" * body
end

"""Decode a 17 character `prefix:body` index string into a HexCell."""
function HexCell(str::AbstractString)
    return HexCell(
        _prefix_to_system(str[1:8]),
        HexFlatTop(_hash_to_hex_coords(str[10:17]))
    )
end

###
# Conversions lon-lat ⇔ index
# (by way of HexCells)
###

"""Compute the lon-lat center of the HexCell defined by an index."""
center(ind::AbstractString) = center(HexCell(ind))

"""Compute the lon-lat vertices of the HexCell defined by an index."""
vertices(ind::AbstractString) = vertices(HexCell(ind))

"""Compute the prefix:body index for the HexCell containing an input lon-lat (in degrees)
in a given HexSystem."""
index(lon, lat, hs::HexSystem) = index(HexCell(lon, lat, hs))

###
# GeoInterface utilities
###

function Polygon(hc::HexCell)
    verts = vertices(hc)[[1:6; 1]] # double up first point to get 6 sides
    return [Point(v...) for v in verts] |> Polygon
end

end # module
