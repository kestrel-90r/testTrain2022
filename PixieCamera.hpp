//-----------------------------------------------
//
//	This file is part of the Siv3D Engine.
//
//	Copyright (c) 2008-2021 Ryo Suzuki
//	Copyright (c) 2016-2021 OpenSiv3D Project
//
//	Licensed under the MIT License.
//
//-----------------------------------------------

# pragma once

namespace s3d
{
	class alignas(16) PixieCamera //:public BasicCamera3D
	{
		//	using BasicCamera3D::BasicCamera3D;

	public:
		// from BasicCamera3D
		static constexpr double DefaultVerticalFOV = 30_deg;
		static constexpr double DefaultNearClip = 0.2;

		virtual ~PixieCamera() = default;
		PixieCamera() = default;
		PixieCamera(const PixieCamera&) = default;
		explicit PixieCamera( const Size& sceneSize, double verticalFOV = DefaultVerticalFOV,
							  const Float3& eyePos = Float3{ 0, 4, -4 },
							  const Float3& focusPos = Float3{ 0, 0, 0 },
							  double nearClip = DefaultNearClip,
							  float cameraRoll = 0.0f,
							  const Float3& orgPos = Float3{ 0, 0, 0 }) noexcept

			: m_sceneSize{ sceneSize }
			, m_verticalFOV{ verticalFOV }
			, m_nearClip{ nearClip }
			, m_eyeLocal{ eyePos }
			, m_focusLocal{ focusPos }
			, m_upDir{ 0, 1, 0 }
			, m_orgLocal{ eyePos }
		{
			setRoll(cameraRoll);
			updateProj();
			updateView();
			updateViewProj();
		}

		void setSceneSize(const Size& sceneSize) noexcept
		{
			m_sceneSize = sceneSize;
			updateProj();
			updateViewProj();
		}

		inline const Size& getSceneSize() const noexcept
		{
			return m_sceneSize;
		}

		void setProjection(const Size& sceneSize, double verticalFOV, double nearClip = DefaultNearClip) noexcept
		{
			m_sceneSize = sceneSize;
			m_verticalFOV = verticalFOV;
			m_nearClip = nearClip;

			updateProj();
			updateViewProj();
		}

		inline const Mat4x4& SIV3D_VECTOR_CALL getProj() const noexcept
		{
			return m_proj;
		}

		void setView(const Float3& eyePos, const Float3& focusPos,
					  const Float3& orgPos = Float3{ 0, 0, 0 },
					  const Float3& upDir = Float3{ 0, 1, 0 }) noexcept
		{
			m_eyeLocal = eyePos;
			m_focusLocal = focusPos;
			m_upDir = upDir;
			m_orgLocal = orgPos;

			updateView();
			updateViewProj();
		}
		inline const Mat4x4& SIV3D_VECTOR_CALL getView() const noexcept
		{
			return m_view;
		}

		PixieCamera& setUpDir(const Float3& upDir) noexcept
		{
			m_upDir = upDir;
			return *this;
		}

		inline const Float3& getUpDir() const noexcept
		{
			return m_upDir;
		}

		inline double getVerticlaFOV() const noexcept
		{
			return m_verticalFOV;
		}

		inline double getNearClip() const noexcept
		{
			return m_nearClip;
		}

		inline const Float3& getOrgEyePos() const noexcept
		{
			return m_orgLocal + m_eyeLocal;
		}
		inline const Float3& getOrgFocusPos() const noexcept
		{
			return m_orgLocal + m_focusLocal;
		}

		inline const Float3& getEyePos() const noexcept
		{
			return m_eyeLocal;
		}

		inline const Float3& getOrgPos() const noexcept
		{
			return m_orgLocal;
		}
		inline const Float3& getFocusPos() const noexcept
		{
			return m_focusLocal;
		}

		Float3 worldToScreenPoint(const Float3& pos) const noexcept
		{
			Float3 v = SIMD_Float4{ DirectX::XMVector3TransformCoord(
									SIMD_Float4{ pos, 0.0f }, m_viewProj) }.xyz();
			v.x += 1.0f;
			v.y += 1.0f;
			v.x *= 0.5f * m_sceneSize.x;
			v.y *= 0.5f;
			v.y = 1.0f - v.y;
			v.y *= m_sceneSize.y;
			return v;
		}

		Float3 screenToWorldPoint(const Float2& pos, const float depth) const noexcept
		{
			Float3 v{ pos, depth };
			v.x /= (m_sceneSize.x * 0.5f);
			v.y /= (m_sceneSize.y * 0.5f);
			v.x -= 1.0f;
			v.y -= 1.0f;
			v.y *= -1.0f;

			const SIMD_Float4 worldPos = DirectX::XMVector3TransformCoord(
										 SIMD_Float4{ v, 0.0f }, m_invViewProj);
			return worldPos.xyz();
		}

		Ray screenToRay(const Vec2& pos) const noexcept
		{
			Float3 orgeye = m_orgLocal + m_eyeLocal;
			const Float3 rayEnd = screenToWorldPoint(pos, 0.1f);
			return Ray{ orgeye, (rayEnd - orgeye).normalized() };
		}

		inline Float3 getLookAtVector() const noexcept
		{
			return (m_focusLocal - m_eyeLocal).normalized();
		}

		Quaternion getLookAtOrientation() const noexcept
		{
			return Quaternion::FromUnitVectorPairs({ Float3::Forward(), Float3::Up() },
													{ getLookAtVector(), Float3::Up() });
		}

		inline const Mat4x4& SIV3D_VECTOR_CALL getInvView() const noexcept
		{
			return m_invView;
		}

		inline const Mat4x4& SIV3D_VECTOR_CALL getViewProj() const noexcept
		{
			return m_viewProj;
		}

		const Mat4x4& SIV3D_VECTOR_CALL getInvViewProj() const noexcept
		{
			return m_invViewProj;
		}

		inline Mat4x4 billboard(const Float3 pos, const Float2 scale) const noexcept
		{
			Mat4x4 m = m_invView;
			m.value.r[0] = DirectX::XMVectorScale(m.value.r[0], scale.x);
			m.value.r[1] = DirectX::XMVectorScale(m.value.r[1], 1.0f);
			m.value.r[2] = DirectX::XMVectorScale(m.value.r[2], scale.y);
			m.value.r[3] = DirectX::XMVectorSet(pos.x, pos.y, pos.z, 1.0f);
			return m;
		}

		void updateProj() noexcept
		{
			const double g = (1.0 / std::tan(m_verticalFOV * 0.5));
			const double s = (static_cast<double>(m_sceneSize.x) / m_sceneSize.y);
			constexpr float e = 0.000001f;

			m_proj = Mat4x4{
				static_cast<float>(g / s), 0.0f, 0.0f, 0.0f,
				0.0f, static_cast<float>(g), 0.0f, 0.0f,
				0.0f, 0.0f, e, 1.0f,
				0.0f, 0.0f, static_cast<float>(m_nearClip * (1.0 - e)), 0.0f
			};
		}

		Quaternion m_eyeRollQ = Quaternion::Identity();
		PixieCamera& setRoll(float cameraRoll) noexcept
		{
			m_eyeRollQ = Quaternion::RotationAxis( m_focusLocal - m_eyeLocal, ToRadians(cameraRoll));
			return *this;
		}

		//カメラは列車座標系で変換されて最終ワールド座標へ
		PixieCamera& updateViewWorld(const Float3& orgpos, const Quaternion& orgrotq) noexcept
		{
			m_upDir = Float3{ 0,1,0 } * m_eyeRollQ;

			Mat4x4 mat = Mat4x4(orgrotq) * Mat4x4::Identity().Translate(orgpos);
			m_eyeWorld = { mat.transformPoint(m_eyeLocal), 0.0f };
			m_focusWorld = { mat.transformPoint(m_focusLocal), 0.0f };

			const SIMD_Float4 upDir{ m_upDir, 0.0f };
			if (m_eyeWorld != m_focusWorld)
				m_view = DirectX::XMMatrixLookAtLH(m_eyeWorld, m_focusWorld, upDir);
			m_invView = m_view.inverse();
			return *this;
		}

		const SIMD_Float4& getEyeWorld() { return m_eyeWorld; }
		const SIMD_Float4& getFocusWorld() { return m_focusWorld; }

		PixieCamera &updateView() noexcept
		{
			m_upDir = Float3{ 0,1,0 } *m_eyeRollQ;

			const SIMD_Float4 eyePos{ (m_orgLocal + m_eyeLocal), 0.0f };
			SIMD_Float4 focusPos{ (m_orgLocal + m_focusLocal), 0.0f };
			const SIMD_Float4 upDir{ m_upDir, 0.0f };
			if ( eyePos != focusPos )
				m_view = DirectX::XMMatrixLookAtLH(eyePos, focusPos, upDir);
			m_invView = m_view.inverse();
			return *this;
		}

		PixieCamera& updateViewProj() noexcept
		{
			m_upDir = Float3{ 0,1,0 } *m_eyeRollQ;
			const SIMD_Float4 upDir{ m_upDir, 0.0f };
			if(m_eyeWorld != m_focusWorld)
				m_view = DirectX::XMMatrixLookAtLH(m_eyeWorld, m_focusWorld, upDir);

			m_viewProj = (m_view * m_proj);
			m_invViewProj = m_viewProj.inverse();
			return *this;
		}

		Mat4x4 m_proj = Mat4x4::Identity();
		Mat4x4 m_view = Mat4x4::Identity();
		Mat4x4 m_invView = Mat4x4::Identity();
		Mat4x4 m_viewProj = Mat4x4::Identity();
		Mat4x4 m_invViewProj = Mat4x4::Identity();

		Size m_sceneSize = Scene::DefaultSceneSize;
		double m_verticalFOV = DefaultVerticalFOV;
		double m_nearClip = DefaultNearClip;

		Float3 m_orgLocal = Float3{ 0, 0, 0 };
		Quaternion m_orgRotQ = Quaternion::Identity();	//orgRotQ x eyePos→ローカル座標

		Float3 m_eyeLocal = Float3{ 0, 4, -4 };			//従来のワールド座標系、視点/注視点
		Float3 m_focusLocal = Float3{ 0, 0, 0 };			//キャラ向き

		SIMD_Float4 m_eyeWorld;
		SIMD_Float4 m_focusWorld;
		Float3 m_upDir = Float3{ 0, 1, 0 };

		Float3 m_eyePosDelta = Float3{ 0, 0, 0 };
		Float3 m_focusPosDelta = Float3{ 0, 0, 0 };

		//列車位置+キャラ位置
		//列車向き+キャラ向き

		// itakawa added camerawork features
		PixieCamera& setOrgPos(const Float3& orgpos) noexcept
		{
			m_orgLocal = orgpos;
			return *this;
		}

		PixieCamera& setOrgRotQ(const Quaternion& orgrot) noexcept
		{
			m_orgRotQ = orgrot;
			return *this;
		}
		Quaternion& getOrgRotQ() noexcept
		{
			return m_orgRotQ ;
		}

		PixieCamera& setEyePos(const Float3& eyepos) noexcept
		{
			m_eyeLocal = eyepos;
			return *this;
		}


		PixieCamera& setEyeX(const float& x) noexcept
		{
			m_eyeLocal.x = x; return *this;
		}
		PixieCamera& setFocusX(const float& x) noexcept
		{
			m_focusLocal.x = x; return *this;
		}
		PixieCamera& setEyeY(const float& y) noexcept
		{
			m_eyeLocal.y = y; return *this;
		}
		PixieCamera& setFocusY(const float& y) noexcept
		{
			m_focusLocal.y = y; return *this;
		}
		PixieCamera& setEyeZ(const float& z) noexcept
		{
			m_eyeLocal.z = z; return *this;
		}
		PixieCamera& setFocusZ(const float& z) noexcept
		{
			m_focusLocal.z = z; return *this;
		}

		PixieCamera& setFocusPos(const Float3& focuspos) noexcept
		{
			m_focusLocal = focuspos;
			return *this;
		}

		inline const Float3& getEyePosDelta() const noexcept
		{
			return m_eyePosDelta;
		}

		inline const Float3& getFocusPosDelta() const noexcept
		{
			return m_focusPosDelta;
		}

		float getBasicSpeed()
		{
			Float3 distance = m_focusLocal - (m_eyeLocal + m_orgLocal);
			Float3 identity = distance.normalized();
			return distance.length() / identity.length();
		}

		//dst/srcの指定で任意座標系を使用
		bool dolly(float forwardback, bool movefocus, bool mask_y = false)
		{
			return dolly(forwardback, movefocus, m_focusLocal, m_eyeLocal, mask_y);
		}

		bool dolly(float forwardback, bool movefocus, const Float3 &dst,const Float3 &src, bool mask_y = false)
		{
			if (forwardback == 0) return false;

			dolly(forwardback, dst, src, mask_y);
			if (movefocus) m_focusLocal += m_focusPosDelta;

			//目標点に達した
//			Float3 dirA = (m_focusLocal - m_eyeLocal);
//			Float3 dirB = (dst - src);
//			return (Math::Sign(dirA.x) != Math::Sign(dirB.x) &&
//					Math::Sign(dirA.y) != Math::Sign(dirB.y) &&
//					Math::Sign(dirA.z) != Math::Sign(dirB.z));
			return false;
		}

		void dolly(float forwardback)
		{
			dolly(forwardback, m_focusLocal, m_eyeLocal);
		}

		//dst/srcの指定で任意座標系を使用
		void dolly(float forwardback, const Float3 &dst, const Float3 &src, bool mask_y = false)
		{
			if (forwardback == 0) return ;
			Float3 dir = (dst - src).normalized();
			m_eyePosDelta = forwardback * dir;

			if (mask_y) m_eyePosDelta.y = 0;

			m_focusPosDelta = m_eyePosDelta;

			m_eyeLocal += m_eyePosDelta;
		}

		void dolly(float forwardback, const Quaternion& rotq)
		{
			if (forwardback == 0) return ;

			SIMD_Float4 axis = rotq.toFloat4();
			axis.setY(0);
			Quaternion q{axis};
			Float3 dir = Mat4x4(q).transformPoint({ 0 , 0,forwardback });

			m_eyeLocal += dir;
			m_focusLocal += dir;
		}

		PixieCamera& panX(float leftright)
		{
			if (leftright == 0) return *this;
			Float3 forward = m_focusLocal - m_eyeLocal;
			Float3 dir = forward.normalized();
			Float3 right = dir.cross(m_upDir);
			m_upDir = right.cross(dir).normalized();

			SIMD_Float4 qrot = DirectX::XMQuaternionRotationAxis(SIMD_Float4{ m_upDir, 0 }, leftright);
			SIMD_Float4 focus = DirectX::XMVector3Rotate(SIMD_Float4{ forward, 0 }, qrot);

			Float3 oldfocus = m_focusLocal;
			m_eyePosDelta = Float3{ 0,0,0 };
			m_focusLocal = m_eyeLocal + focus.xyz();
			m_focusPosDelta = m_focusLocal - oldfocus;
			return *this;
		}

		PixieCamera& panY(float updown)
		{
			if (updown == 0) return *this;
			Float3 forward = m_focusLocal - m_eyeLocal;
			Float3 dir = forward.normalized();
			Float3 right = dir.cross(m_upDir);
			m_upDir = right.cross(forward).normalized();

			SIMD_Float4 qrot = DirectX::XMQuaternionRotationAxis(SIMD_Float4{ right, 0 }, updown);
			SIMD_Float4 focus = DirectX::XMVector3Rotate(SIMD_Float4{ forward, 0 }, qrot);

			Float3 oldfocus = m_focusLocal;
			m_eyePosDelta = Float3{ 0,0,0 };
			m_focusLocal = m_eyeLocal + focus.xyz();
			m_focusPosDelta = m_focusLocal - oldfocus;
			return *this;
		}

		PixieCamera& tilt(float updown)
		{
			panY(updown);
			return *this;
		}

		//dst/srcの指定で任意座標系を使用
		PixieCamera& trackX(float leftright, bool mask_y = false )
		{ return trackX(leftright, m_focusLocal, m_eyeLocal, mask_y); }

		PixieCamera& trackX(float leftright, const Float3 &dst, const Float3 &src, bool mask_y = false)
		{
			if (leftright == 0) return *this;
			Float3 fwd = dst - src;
			Float3 dir = fwd.normalized();
			Float3 right = dir.cross(m_upDir).normalized();

			dir = leftright * right;
			if (mask_y) dir.y = 0;
			m_eyePosDelta = dir;
			m_focusPosDelta = dir;
			m_eyeLocal += dir;
			m_focusLocal += dir;
			return *this;
		}

		//dst/srcの指定で任意座標系を使用
		PixieCamera& craneY(float updown)
		{ return craneY( updown, m_focusLocal, m_eyeLocal); }

		PixieCamera& craneY(float updown, const Float3& dst, const Float3& src)
		{
			if (updown == 0) return *this;
			Float3 fwd = dst - src;
			Float3 dir = fwd.normalized();
			Float3 right = dir.cross(m_upDir);
			m_upDir = right.cross(dir).normalized();

			dir = updown * m_upDir;
			m_eyePosDelta = dir;
			m_focusPosDelta = dir;
			m_eyeLocal += dir;
			m_focusLocal += dir;
			return *this;
		}

		PixieCamera& arcX(float leftright)
		{
			if (leftright == 0) return *this;
			Float3 forward = m_eyeLocal - m_focusLocal;
			Float3 dir = (-forward).normalized();
			Float3 right = m_upDir.cross(dir);
			m_upDir = dir.cross(right).normalized();

			SIMD_Float4 qrot = DirectX::XMQuaternionRotationAxis(SIMD_Float4{ m_upDir, 0 }, leftright);
			SIMD_Float4 eye = DirectX::XMVector3Rotate(SIMD_Float4{ forward, 0 }, qrot);

			Float3 oldeye = m_eyeLocal;
			m_eyeLocal = eye.xyz() + m_focusLocal;
			m_eyePosDelta = m_eyeLocal - oldeye;
			m_focusPosDelta = Float3{ 0,0,0 };
			return *this;
		}

		PixieCamera& arcY(float updown)
		{
			if (updown == 0) return *this;
			Float3 forward = m_eyeLocal - m_focusLocal;
			Float3 dir = (-forward).normalized();
			Float3 right = m_upDir.cross(dir);
			m_upDir = dir.cross(right).normalized();

			SIMD_Float4 qrot = DirectX::XMQuaternionRotationAxis(SIMD_Float4{ right, 0 }, updown);
			SIMD_Float4 eye = DirectX::XMVector3Rotate(SIMD_Float4{ forward, 0 }, qrot);
			Float3 oldeye = m_eyeLocal;
			m_eyeLocal = eye.xyz() + m_focusLocal;
			m_eyePosDelta = m_eyeLocal - oldeye;
			m_focusPosDelta = Float3{ 0,0,0 };
			return *this;
		}

		//カメラ行動制限範囲
//		Box boxCollider = {};
//		PixieCamera& setBoxCollider( Box box ) { boxCollider = box; }
//		Box getBoxCollider( Box box) { return boxCollider; }

		//移動座標系
		Float3 moveForward = {0,0,1};
		Float3 moveUp = { 0,1,0 };
		Float3 moveRight = { 1,0,0 };
		Quaternion moveQRotate = Quaternion::Identity();

		Float3 getF() { return moveForward; }
		Float3 getU() { return moveUp; }
		Float3 getR() { return moveRight; }
		Quaternion getQ() { return moveQRotate; }

		PixieCamera& setMovingCoordinate(const Float3& dst, const Float3& src)
		{
			Float3& forward = moveForward;
			Float3& up = moveUp;
			Float3& right = moveRight;

			if (dst == src)
			{
				forward = { 0,0,1 };
				up = { 0,1,0 };
				right = { 1,0,0 };
				moveQRotate = Quaternion::Identity();
				return *this;
			}

			Float3 _up = Float3{ 0,1,0 };
			_up = up;

			Float3 _forward = dst - src;
			Float3 _dir = _forward.normalized();
			if (_dir == _up) _dir += Float3{ 0.000001,0,0 };	//上と方向が同じ
			Float3 _right = _dir.cross(_up).normalized();
			_up = _dir.cross(_right).normalized();

			forward = _dir;
			up = _up;
			right = _right;
			
			moveQRotate = Quaternion::FromUnitVectorPairs(
				{ Float3{0,0,1}, Float3{0,1,0} }, { -_dir, -_up });
			return *this;
		}

		static Mat4x4 getLookAt(const Float3& dst, const Float3& src)
		{
			if (dst == src) return Mat4x4::Identity();
			Float3 _up = Float3{ 0,1,0 };

			Float3 _front = (dst - src).normalized();
			if (_front == _up) _front += Float3{ 0.000001,0,0 };	//上と方向が同じ

			Float3 _right = _front.cross(_up).normalized();
			_up = _front.cross(_right).normalized();
			Quaternion _q = Quaternion::FromUnitVectorPairs({ Float3{0,0,1}, Float3{0,1,0} }, { -_front, -_up });

			Mat4x4 mat = Mat4x4::Set(
				_front.x, _front.y, _front.z, _right.x,
				_front.x, _front.y, _front.z, _right.y,
				_front.x, _front.y, _front.z, _right.z,
				_up.x, _up.y, _up.z, _q.getW() );

			return mat;
		}


		static Quaternion getQLookAt(const Float3& dst, const Float3& src,
							   Float3* up = nullptr, Float3* right = nullptr, Float3* foward = nullptr)
		//ポインタを指定する場合は結果を引数で取得
		{
			if (dst == src) return Quaternion::Identity();
			Float3 _up = Float3{ 0,1,0 };
			if (up != nullptr) _up = *up;

			Float3 _forward = dst - src;
			Float3 _dir = _forward.normalized();
			if (_dir == _up) _dir += Float3{ 0.000001,0,0 };	//上と方向が同じ

			Float3 _right = _dir.cross(_up).normalized();
			_up = _dir.cross(_right).normalized();

			if (up != nullptr)    *up = _up;
			if (right != nullptr) *right = _right;
			if (foward != nullptr) *foward = _dir;

			return Quaternion::FromUnitVectorPairs(
				{ Float3{0,0,1}, Float3{0,1,0} },{ -_dir, -_up });
		}

		Quaternion &setQLookAtOrg(const Float3& dst, const Float3& src)
			//ポインタを指定する場合は結果を引数で取得
		{
			if (dst == src) return m_orgRotQ = Quaternion::Identity();
			Float3 _up = Float3{ 0,1,0 };

			Float3 _forward = dst - src;
			Float3 _dir = _forward.normalized();
			if (_dir == _up) _dir += Float3{ 0.000001,0,0 };	//上と方向が同じ

			Float3 _right = _dir.cross(_up).normalized();
			_up = _dir.cross(_right).normalized();

			return m_orgRotQ = Quaternion::FromUnitVectorPairs(
				{ Float3{0,0,1}, Float3{0,1,0} },{ -_dir, -_up });
		}

/*
		// dstはワールド座標指定
		Quaternion getQLookAtWorld(const Float3& dst,
							  Float3* up = nullptr, Float3* right = nullptr) const
		{
			if (dst == (m_eyePos + m_orgPos)) return Quaternion::Identity();
			Float3 _up = Float3{ 0,1,0 };
			if (up != nullptr) _up = *up;

			Float3 _forward = dst - (m_eyePos + m_orgPos);
			Float3 _dir = _forward.normalized();
			if (_dir == _up) _dir += Float3{ 0.000001,0,0 };	//上と方向が同じ
			Float3 _right = _dir.cross(_up).normalized();
			_up = _dir.cross(_right).normalized();
			if (up != nullptr)    *up = _up;
			if (right != nullptr) *right = _right;

			return Quaternion::FromUnitVectorPairs({ Float3{0,0,1}, Float3{0,1,0} },
													{ -_dir, -_up });
		}
		// dstはローカル座標指定m_orgPos原点
		Quaternion getQLookAtLocal(const Float3& dst,
							  Float3* up = nullptr, Float3* right = nullptr) const
		{
			if (dst == m_eyePos) return Quaternion::Identity();
			Float3 _up = Float3{ 0,1,0 };
			if (up != nullptr) _up = *up;

			Float3 _forward = dst - m_eyePos;
			Float3 _dir = _forward.normalized();
			if (_dir == _up) _dir += Float3{ 0.000001,0,0 };	//上と方向が同じ
			Float3 _right = _dir.cross(_up).normalized();
			_up = _dir.cross(_right).normalized();
			if (up != nullptr)    *up = _up;
			if (right != nullptr) *right = _right;

			return Quaternion::FromUnitVectorPairs({ Float3{0,0,1}, Float3{0,1,0} },
													{ -_dir, -_up });
		}
*/

		Mat4x4 getRotMat() const
		{
			return Mat4x4{ getQLookAt(m_focusLocal, m_eyeLocal) };
		}

		Quaternion getQForward() const
		{
			return getQLookAt(m_focusLocal, m_eyeLocal);
		}

		Quaternion getQUp() const
		{
			return getQLookAt(m_upDir, m_eyeLocal);
		}

		Quaternion getQRight() const
		{
			Float3 forward = (m_focusLocal - m_eyeLocal).normalized();
			Float3 right = forward.cross(m_upDir).normalized();
			return getQLookAt(right, m_eyeLocal);
		}

		// posはワールド座標指定
		Mat4x4 getMatInfront(Float3 pos, Float2 scale = Float2{ 1,1 }, float distance = 5.0f) const
		{
			Float3 infront = (m_eyeLocal + m_orgLocal - pos).normalized();
			pos = infront * distance;
			return billboard(pos, scale);
		}

	};
}
