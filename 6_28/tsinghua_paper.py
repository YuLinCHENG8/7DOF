



def _solve_theta4_from_triangle(S: np.ndarray, W: np.ndarray, p: NeroParams) -> Optional[float]:
    l_sw = np.linalg.norm(W - S)
    l_se = abs( p.d_i[ 2])
    l_ew = abs( p.d_i[ 4])
    c4 = (l_sw**2 - l_se**2 - l_ew**2) / (2.0 * l_se * l_ew)
    c4 = np.clip(c4, -1.0, 1.0)
    return  math.acos(c4)

def _compute_swe_from_target(T07: np.ndarray, p: NeroParams) -> Tuple[ np.ndarray, np.ndarray,  Optional[float], np.ndarray]:
    R = T07[:3, :3]
    p_target = T07[:3, 3]
    z7 = R[:, 2]
    d6 = float( p.d_i[ 6])
    d1 = float( p.d_i[ 0])

    # 末端法兰中心
    O7 = p_target - p.post_transform_d8 * z7
    # 腕心 W：从法兰再回退 d6
    W = O7 - d6 * z7
    # 肩心 S：固定在基座上方 d1
    S = np.array([ 0.0, 0.0, d1], dtype=float)

    # 余弦定理求 θ4 绝对值
    q4_abs = _solve_theta4_from_triangle(S, W, p)

    # 肩→腕单位向量
    v_sw = W - S
    n_sw = np.linalg.norm(v_sw)
    u_sw = v_sw / n_sw if n_sw > 1e-12else  np.array([ 0.0, 0.0, 1.0])

    return S, W, q4_abs, u_sw



def _elbow_from_arm_angle(S: np.ndarray, W: np.ndarray, theta0: float, p: NeroParams) -> Optional[ np.ndarray]:
    l_se = abs( p.d_i[ 2])
    l_ew = abs( p.d_i[ 4])
    sw = W - S
    l_sw = np.linalg.norm(sw)
    u_sw = sw / l_sw

    # 圆心 C 在 SW 线上的投影
    x = (l_se**2 - l_ew**2 + l_sw**2) / (2.0 * l_sw)
    r2 = l_se**2 - x**2
    r = math.sqrt( max(0.0, r2))
    C = S + x * u_sw

    # 构造圆平面坐标系 e1, e2
    os_vec = S.copy()
    t = np.cross(os_vec, u_sw)
    e1 = t / np.linalg.norm(t)
    e2 = np.cross(u_sw, e1)
    e2 = e2 / np.linalg.norm(e2)

    # 由臂角 theta0 计算 E
    E = C + r * ( math.cos(theta0) * e1 + math.sin(theta0) * e2)
    return E




def _solve_q123_from_swe(E: np.ndarray, W: np.ndarray, q4: float, p: NeroParams) -> List[ np.ndarray]:
    d0 = p.d_i[ 0]
    d2 = p.d_i[ 2]
    d4 = p.d_i[ 4]
    Ex, Ey, Ez = E

    # q2
    c2 = (Ez - d0) / d2
    c2 = np.clip(c2, -1.0, 1.0)
    s2_abs = math.sqrt( max(0.0, 1.0 - c2**2))

    s4 = math.sin(q4)
    c4 = math.cos(q4)
    sols = []

    # 遍历 s2 正负两种构型
    for s2 in (s2_abs, -s2_abs):
        # q1
        c1 = -Ex / (d2 * s2)
        s1 = -Ey / (d2 * s2)
        n1 = math.hypot(c1, s1)
        c1 /= n1
        s1 /= n1
        q1 = math.atan2(s1, c1)
        q2 = math.atan2(s2, c2)

        # q3
        v = W - E
        col2 = -v / d4
        u1, u2, u3 = col2
        b1 = (s2 * c1 * c4 - u1) / s4
        b2 = (u2 - s1 * s2 * c4) / s4
        s3 = s1 * b1 + c1 * b2
        c2c3 = -c1 * b1 + s1 * b2
        c3 = c2c3 / c2 ifabs(c2) > 1e-8else (u3 + c2 * c4) / (s2 * s4)
        n3 = math.hypot(s3, c3)
        s3 /= n3
        c3 /= n3
        q3 = math.atan2(s3, c3)

        sols.append(np.array([q1, q2, q3]))
    return sols



def _extract_567_from_T47_paper(T47: np.ndarray ) -> List[ np.ndarray]:
    sols = []
    c6 = np.clip(T47[ 1, 2], -1.0, 1.0)
    for sgn in (1.0, -1.0):
        s6 = sgn * math.sqrt( max(0.0, 1.0 - c6**2))
        if abs(s6) < 1e-8:
            continue
        th6 = math.atan2(s6, c6)
        th5 = math.atan2(T47[ 2, 2] / s6, T47[0, 2] / s6)
        th7 = math.atan2(T47[ 1, 1] / s6, -T47[1, 0] / s6)
        sols.append(np.array([th5, th6, th7]))
    return sols

def _get_theta0_feasible_region(T07: np.ndarray, p: NeroParams, step: float = 0.01) -> List[float]:
    feasible = []
    for theta0 in  np.arange(-math.pi, math.pi, step):
        if _ik_one_arm_angle(T07, theta0, p):
            feasible.append( float(theta0))
    return feasible


def _weight_limits(q: float, q_min: float, q_max: float) -> float:
    span = q_max - q_min
    x = 2.0 * (q - (q_min + q_max)*0.5) / span
    a = 2.38
    b = 2.28
    if x >= 0:
        den = math.exp(a*( 1-x)) - 1
        return b*x / den
    else:
        den = math.exp(a*( 1+x)) - 1
        return -b*x / den
    


def _optimal_theta0(feasible_theta0, T07, p, q_prev):
    best_cost = inf
    best_t = feasible_theta0[0]
    for t in feasible_theta0:
        sols = _ik_one_arm_angle(T07, t, p)
        for q_full in sols:
            q = q_full[:7]
            cost = 0
            for i inrange(7):
                lo, hi = p.joint_limits[i]
                w = _weight_limits(q[i], lo, hi)
                dq = abs(q[i] - q_prev[i])
                cost += w * dq*dq
            if cost < best_cost:
                best_cost = cost
                best_t = t
    return best_t

    