{-|
Module      : OrbitalMechanics
Description : This module contains basic orbital mechanic functions.
License     : CC0
Maintainer  : frosch03@frosch03.de
Stability   : experimental
-}
module System.KSP.OrbitalMechanics
where

import System.KSP.DataConstructors
import System.KSP.Datatype.Body (DeltaV)

type Gravity = Double
type Isp = Double

-- | 'Radius' of 'Orbit' at the current Position
type Radius = Double 

-- | 'Speed' is an alias for 'Double'
type Speed  = Double

type Mass  = Double

-- | 'var_G' [Nm^2/kg^2] is the Gravitation constant in newton meter
-- squared over kilo gramms squared
var_G :: GravConst
var_G = 6.674e-11

-- | 'g_ker' [m/s^2] is acceleration due to gravity near Kerbin surface
g_ker = 9.80665

-- | 'semiMajor' calculates the semi major axis of an orbit.
semiMajor   :: Orbit Body -> Double
semiMajor o
    = 0.5 * ((apoapsis o) + (periapsis o) + (2 * (r . celestial . centerBody $ o)))


-- | 'v' takes an orbit and a radius (from the center of the
-- centerBody) and calculates the orbital speed at that position.
v :: Orbit Body -> Radius -> Speed
v o radius =  sqrt $ mue * ( (2/radius) - (1/(semiMajor o)) ) 
    where
      mue = var_G * (m . celestial . centerBody $ o)


-- | 'v_e' calculates the escape velocity of that body.
v_e :: Body -> Speed
v_e b =  sqrt $ 2 * mue / (r . celestial $ b) 
    where
      mue = var_G * (m . celestial $ b)

      

-- | 'hohmann' takes two orbits around the same centerBody. It
-- calculates both (v1 and v2 ) delta V changes for a hohmann
-- transfair.
hohmann :: Orbit Body -> Orbit Body -> (Double, Double)
hohmann o1 o2
    |   centerBody o1 /= centerBody o2
      || o1 == o2
    = (0,0)
hohmann o1 o2
    = ( sqrt(mue/r1) * (sqrt((2*r2)/(r1+r2)) - 1)
      , sqrt(mue/r2) * (1 - (sqrt((2*r1)/(r1+r2))))
      )
    where
      mue = var_G * (m . celestial . centerBody $ o1)
      r1 = (r . celestial . centerBody $ o1) + (apoapsis o1)
      r2 = (r . celestial . centerBody $ o2) + (apoapsis o2)


-- | 'tsky' computes Tsiolkovsky delta v given specific impulse I_sp,
-- standard gravity g as well as initial and final masses m_0 and m_f.
tsky :: Gravity -> Isp -> Mass -> Mass -> DeltaV
tsky g isp m_0 m_f = (g*isp)*log(m_0 / m_f)

-- | `tskyKer` is a convenience function for 'tsky' with
-- Kerbin acceleration assumed
tskyKer :: Isp -> Mass -> Mass -> DeltaV
tskyKer isp m_0 m_f = (g_ker*isp)*log(m_0 / m_f)
