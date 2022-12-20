# include <Siv3D.hpp> // OpenSiv3D v0.6.5

typedef uint32 uint;

# include "PixieMesh.hpp"
# include "PixieCamera.hpp"
# include "LineString3D.hpp"
# include "Maple2.hpp"

constexpr ColorF        GRAY1 = { 0.01, 0.01, 0.01, 0.5 };
constexpr ColorF        GRAY3 = { 0.03, 0.03, 0.03, 0 };
constexpr ColorF        GRAY5 = { 0.05, 0.05, 0.05, 0 };
constexpr ColorF        GRAY8 = { 0.08, 0.08, 0.08, 0 };

constexpr ColorF        BGCOLOR = GRAY5;
constexpr TextureFormat TEXFORMAT = TextureFormat::R8G8B8A8_Unorm_SRGB;
constexpr Size			WINDOWSIZE = { 1440, 768 };

constexpr StringView	APATH = U"../Asset/";


//							blend.opAlpha = BlendOp::Max;\

#define TRANSPARENTBG(RTEX)	BlendState blend = BlendState::Default2D;\
							blend.srcAlpha = Blend::SrcAlpha;\
							blend.dstAlpha = Blend::DestAlpha;\
							blend.opAlpha = BlendOp::Max;\
							const ScopedRenderStates3D bs{ BlendState::OpaqueAlphaToCoverage,\
								  SamplerState::RepeatAniso, RasterizerState::SolidCullBack };\
							const ScopedRenderTarget3D rtsub{ (RTEX).clear(ColorF{ GRAY1 }) }

#define RENDERTEXTURE(RTEX,POS)	\
		Graphics3D::Flush();\
		(RTEX).resolve();\
		Shader::LinearToScreen((RTEX), (POS))


template <typename V, typename T>
inline V guiValue(const Font& font, const String& text, Rect rect, float wheel, V initvalue, V value, V& prev, T min, T max)
{
	if (rect.mouseOver() && wheel != 0) value = Math::Clamp(value + wheel, min, max);
	if (rect.leftClicked()) value = initvalue;
	font(text).draw(rect, rect.mouseOver() ? ColorF(Palette::Red) : ColorF(Palette::White));
	if (prev != value) { LOG_INFO(text); prev = value; }
	return value;
}
template <typename V>
inline V guiValue(const Font& font, const String& text, Rect rect, V value, V& prev)
{
	font(text).draw(rect, rect.mouseOver() ? ColorF(Palette::Red) : ColorF(Palette::White));
	if (prev != value) { LOG_INFO(text); prev = value; }
	return value;
}

template <typename V>
inline V AddLooped(V value, const V& max,const V& amount)
{
	value += amount;
	if (value > max) value -= max;
	return value;
}

const Font* fontS;
const Font* fontXS;
PixieCamera cameraMain(WINDOWSIZE);
PixieCamera cameraSub(WINDOWSIZE);

enum { CAR0, CAR1, CAR2, CAR3, CAR4, CAR5, NUM_CARS };
enum { NPC0, NPC1, NPC2, USER, NPC4, NPC5, NPC6, NUM_TRAINS };

struct Train
{
	double progressT = 0.25;
	double spacing = 0.00394;
	double speedT = 0.000001;

	Array<Float3> carPos;
	Array<Quaternion> carRotQ;

	Train()
	{
		carPos.resize(NUM_CARS, { 0,0,0 });
		carRotQ.resize(NUM_CARS, Quaternion::Identity());
	}
};


Array<Train> Trains(NUM_TRAINS);
PixieMesh meshMapleSub;
int atCar = CAR1;


//タグ定義
enum {
	TRAIN_TAIL = 29, TRAIN_MID = 28, TRAIN_HEAD = 27,
	COLLIDER_TAIL = 32, COLLIDER_MID = 31, COLLIDER_HEAD = 30,
	ARROW = 33, CAMERA = 34,
};

Array<Entity> entMaple = {
	//地下鉄駅
	Entity{ "A1",  0}, Entity{ "D6", 1 }, Entity{ "A2", 2 }, Entity{ "A3", 3 },
	Entity{ "B1",  4}, Entity{ "B2", 5 }, Entity{ "B3", 6 }, Entity{ "C1", 7 },
	Entity{ "C2",  8}, Entity{ "C3", 9 }, Entity{ "C4", 10}, Entity{ "C5", 11},
	Entity{ "C6", 12}, Entity{ "D1", 13}, Entity{ "D2", 14}, Entity{ "D3", 19},
	Entity{ "D4", 15}, Entity{ "D5", 16}, Entity{ "ES", 17},

	//地下鉄線路			  地下鉄列車	
	Entity{ "SW", 18 },
	Entity{ "R1", 20 }, Entity{ "R2", 21}, Entity{ "R3", 22},
	Entity{ "S1", 23 }, Entity{ "S2", 24}, Entity{ "S3", 25},

	//ドローン		   列車(後尾車両)  	           列車(車両)		          列車(先頭車両)
	Entity{ "DR", 26}, Entity{ "TR1", TRAIN_HEAD },Entity{ "TR2", TRAIN_MID },Entity{ "TR3", TRAIN_TAIL },

	//コライダー   後尾車両        先頭                          車両		            
	Entity{ "CL1", COLLIDER_TAIL}, Entity{ "CL2", COLLIDER_MID}, Entity{ "CL3", COLLIDER_HEAD},

	//矢印              カメラ
	Entity{ "AR", ARROW},	Entity{ "CAM", CAMERA},
};

void updateMainCamera(PixieCamera& camera, const Float3& orgpos, const Quaternion& orgrotq,
					   const Maple& maple, const OrientedBox& ob)
{
	Float2 delta = Cursor::DeltaF();
	Float2 point2D = Cursor::PosF();
	Float3 distance = {};

	float speedM = 1;
	if (MouseM.down())								//中ボタンドラッグ：回転
	{
		const Ray mouseRay = camera.screenToRay(Cursor::PosF());
		if (mouseRay.direction.hasNaN())
			assert(1);
	}

	if (MouseM.pressed())
	{
		camera.panX(delta.x * speedM / 100);
		camera.panY(delta.y * -speedM / 100);
	}

	if (MouseR.pressed())							//右ボタンドラッグ：平行移動
	{
		camera.trackX(delta.x * speedM / 20);
		camera.craneY(delta.y * speedM / 20);
	}

	if (KeyLControl.pressed() || MouseM.pressed())	//中ボタンウィール：拡縮
	{
		Float3 eyepos = cameraMain.getOrgEyePos();
		Float3 focuspos = cameraMain.getFocusPos();
		bool rev = camera.dolly((float)Mouse::Wheel() * speedM / 5, true);	//中ボタンウィール：拡縮
		if (rev)
		{
			rev = camera.setEyePos(eyepos).dolly((float)Mouse::Wheel() * speedM / 100, true);
			if (rev)
			{
				rev = camera.setEyePos(eyepos).dolly((float)Mouse::Wheel() * speedM / 1000, true);
				if (rev) camera.setEyePos(eyepos);
			}
		}
	}

	float speedK = camera.getBasicSpeed() / 10000;
	Float3 oldeyepos = camera.getEyePos();
	Float3 oldfocuspos = camera.getFocusPos();

	if (KeyW.pressed()) camera.dolly(+speedK, true, true);
	if (KeyS.pressed()) camera.dolly(-speedK, true, true);
	if (KeyA.pressed()) camera.trackX(+speedK, true);
	if (KeyD.pressed()) camera.trackX(-speedK, true);

	if ('#' == maple.GetChar(Float2{ camera.getEyePos().x + ob.w / 2, oldeyepos.z }, ob))
		camera.setEyeX(oldeyepos.x).setFocusX(oldfocuspos.x);

	if ('#' == maple.GetChar(Float2{ oldeyepos.x + ob.w / 2, camera.getEyePos().z }, ob))
		camera.setEyeZ(oldeyepos.z).setFocusZ(oldfocuspos.z);

	camera.setUpDir(Float3{ 0,1,0 })
		.updateViewWorld(orgpos, orgrotq)
		.updateViewProj();

}


void testLine(const LineString3D& railway, PixieCamera& camera, PixieMesh& mesh)
{
	Box b{ 0,0,0 , 0.5 };
	for (int i = 0; i < railway.size(); i += 1)
		b.movedBy(railway[i]).drawFrame(Palette::Cyan);

	mesh.setMove(camera.getEyeWorld().xyz() + Float3{ 0,+2,0 })
		.setRotateQ(camera.getQForward() * Trains[USER].carRotQ[atCar])
		.drawMesh(ARROW, Palette::Green);
}


void updateTrain(LineString3D& railway, const double progressT, const double spacing,
				 Array<Float3>& carPos, Array<Quaternion>& carRotQ)
{
	for(int i=0; i<carPos.size(); i++ )
		carPos[i] = railway.getCurvedPoint( AddLooped(progressT,1.0, spacing * (i+2) ) );

	for (int i=0; i < carPos.size()-1; i++)
		carRotQ[i] = PixieCamera::getQLookAt(carPos[i], carPos[i+1]);
}

void changeCar( int car, int meshid, PixieMesh& meshMaple, const Array<Float3>& carPos,
	          const Array<Quaternion>& carRotQ)
{
	// 乗車列車コライダーを取得
	OrientedBox ob = meshMaple.getCollider(meshid);
	ob = Geometry3D::TransformBoundingOrientedBox(ob, meshMaple.getMat());
	if (ob.contains(cameraMain.getEyeWorld().xyz()))
	{
		//列車コライダー変更検知
		if (atCar != car)
		{
			//視点座標を相対で保持
			Float3 dir = cameraMain.getFocusPos() - cameraMain.getEyePos();

			//逆行列でレール回転戻して新列車のローカル座標を算出
			Float3 eyepos = Mat4x4(carRotQ[car]).inverse()
				.transformPoint(cameraMain.getEyeWorld().xyz() - carPos[car]);

			cameraMain.setEyePos(eyepos);
			cameraMain.setFocusPos(eyepos + dir);
			atCar = car;
		}
	}

	if (atCar == car) meshMapleSub = meshMaple;
}


void drawTrain(PixieMesh& meshMaple, const Array<Float3>& carPos, const Array<Quaternion>& carRotQ)
{
	const Array<int> car{ TRAIN_TAIL, TRAIN_MID, TRAIN_MID, TRAIN_MID, TRAIN_HEAD };

	//5両分の列車を描画、乗り換え処理
	for (int i = 0; i < car.size(); i++)
	{
		meshMaple.setMove(carPos[i]).setRotateQ(carRotQ[i]).drawMesh( car[i] );
		changeCar(i, car[i], meshMaple, carPos, carRotQ);
	}

	//列車位置を曲線座標から更新
	cameraMain.setOrgPos( carPos[atCar] );
}

void Main()
{
	const Font FONTS(15); fontS = &FONTS;
	const Font FONTXS(8, FileSystem::GetFolderPath(SpecialFolder::SystemFonts) + U"MSGOTHIC.TTC"); fontXS = &FONTXS;
	static MSRenderTexture rtexMain = { (unsigned)WINDOWSIZE.x, (unsigned)WINDOWSIZE.y, TEXFORMAT, HasDepth::Yes };
	static MSRenderTexture rtexSub = { (unsigned)WINDOWSIZE.x, (unsigned)WINDOWSIZE.y, TEXFORMAT, HasDepth::Yes };

	//ウィンドウ初期化
	Window::Resize(WINDOWSIZE);
	Window::SetStyle(WindowStyle::Sizable);

	const Size SWPOS = Window::GetState().virtualSize - Point{ 805,365 };
	const RectF SUBWINDOW = RectF{ SWPOS, 800,360 };
	const RectF MAINWINDOW = RectF{ 0,0, WINDOWSIZE };

	constexpr double RADIUS = 1000;
	constexpr double VOLATILITY = 500;

	// --- --- --- 

	const Array<std::string> mapSubway1 {
	//  0m  4m  8m  12m 16m 20m 24m 
		"R3..R2..R2..R2..R2..S2..R1..\n"
		":...:...:...:...:...:...:...\n",

		"SW..SW..SW..SW..SW..SW..SW..\n"
		":...:...:...:...:...:...:...\n",

		"R3..S3..R2..R2..R2..R2..R1..\n"
		":...:...:...:...:...:...:...\n",

		"R3..R2..R2..R2..R2..R2..R1..\n"
		":...:...:...:...:...:...:...\n",

		":...R3..R2..R2..R2..R1..:...\n"
		":...:...:...:...:...:...:...\n",

		":...:...R3..R2..R1..:...:...\n"
		":...:...:...:...:...:...:...\n",
	};

	//列車室内のコライダー
	std::string mapInterior1 =
	{//  0m  1m  2m  3m
		"####:####\n"
		"####:####\n"
		"##..:..##\n"
		"##..:..##\n"
		"#...:...#\n"
		"#...:...#\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"#...:...#\n"
		"#...:...#\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"#...:...#\n"
		"#...:...#\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"##..:..##\n"
		"#...:...#\n"
		"#...:...#\n"
		"##..:..##\n"
		"##..:..##\n"
		"####:####\n"
		"####:####\n"
	};

	Maple mapInterior{ mapInterior1 };

	// 環状線計算
	LineString3D trackway;
	for (int r = 0; r < 10; r++)
	{
		constexpr Vec2 CENTER{ 0,0 };
		const Float2 pos2 = CENTER + Circular(RADIUS, ToRadians(r * 36));
		Float3 pos = Float3{ pos2.x + VOLATILITY * (Random() - 0.5),
									  VOLATILITY * Random(),
							 pos2.y + VOLATILITY * (Random() - 0.5) };
		trackway.emplace_back(pos);
	}

	LineString3D curvedway = trackway.catmullRomClosed(80);						//曲線化
	double len = curvedway.updateLengthList();

	LineString3D railway{ curvedway.getEquatePoints((int)(len / 4 )) };			//等速化
	railway.updateLengthList();

	Array<LineString3D> Railways(NUM_TRAINS);									//複線化
	for (int i = 0; i < railway.size(); i++)
	{
		Float3 right;
		Vec3 org = (i == 0) ? railway.back() : railway[i - 1];

		PixieCamera::getQLookAt(org, railway[i], nullptr, &right);
		Railways[0].emplace_back(org - 4 * right * 3);
		Railways[1].emplace_back(org - 4 * right * 2);
		Railways[2].emplace_back(org - 4 * right);
		Railways[3].emplace_back(org);
		Railways[4].emplace_back(org + 4 * right);
		Railways[5].emplace_back(org + 4 * right * 2);
		Railways[6].emplace_back(org + 4 * right * 3);
	}

	for (int i = 0; i < Railways.size(); i++)									//車両間隔
	{
		Railways[i].updateLengthList();

		Float3 p0 = Railways[i].getCurvedPoint(0.5);
		for (int s = 100; s < 200; s++)
		{
			double spacing = double(s) / 40000;
			Float3 p1 = Railways[i].getCurvedPoint( 0.5+spacing );
			double l = p1.distanceFrom(p0);
			if ( l >= 23.9f )
			{
				Trains[i].spacing = spacing;
				break;
			}
		}
	}

	//2kmの地下鉄路線(単線)
	std::string subwaymap = {};

	//最大長は曲線の長さSizeに制限
	int b = 0;
	for (int i = 0; i < int(len+1); i ++)
	{
		if (i % 200 == 0) b = Random(mapSubway1.size()-1);
		subwaymap += mapSubway1[b];
	}

	Maple mapSubway = { subwaymap };
	Chunk chunkSubway;

	constexpr Float4 mainCamEye = { 10, 1, 0, 0 };
	constexpr Float3 mainCamFocus = { 0, 1, 0 };
	cameraMain = PixieCamera(WINDOWSIZE, 45_deg, mainCamEye.xyz(), mainCamFocus, 0.05, mainCamEye.w);

	constexpr Float4 subCamEye = mainCamEye + Float4{ 10,  0, 0, 0 };
	cameraSub = PixieCamera(WINDOWSIZE, 45_deg, subCamEye.xyz(), mainCamEye.xyz(), 0.05, 0.0);

	PixieMesh meshMaple;
	meshMaple.initModel(APATH + U"Subway.07.glb", MODELNOA, WINDOWSIZE, USE_INDIVIDUAL);

	//車両のコライダー調整
	meshMaple.obCollider[TRAIN_HEAD].size += { 0, 0, -1.0};
	meshMaple.obCollider[TRAIN_MID].size += { 0, 0, -1.0};
	meshMaple.obCollider[TRAIN_TAIL].size += { 0, 0, -1.0};

	// 全EntityList生成
	initMaple(mapSubway, entMaple);
	Array<Float3> vec(mapSubway.size.x);

	constexpr int SIZE = 4;
	auto& DL = mapSubway.displaylist;

	for (uint i = 0; i < railway.size(); i++)
	{
		Float3 right = {};

		Vec3 src = (i == 0) ? railway.back() : railway[i-1];
		Vec3 dst = railway[i];

		Quaternion q = PixieCamera::getQLookAt(dst, src, nullptr, &right);

		uint j = i * mapSubway.size.x;
		auto ref = mapSubway.refgrid[j];

		for (uint x = 0; x < mapSubway.size.x; x++)
		{
			ref = mapSubway.refgrid[j + x];
			if (ref != NOENTITY)
			{
				DL[ref].qrotate = q;
				DL[ref].pos = src + right * DL[ref].pos.x * SIZE;
			}
		}
	}

	cameraMain.setEyePos({ 0,1,5 });
	cameraMain.setFocusPos(cameraMain.getEyePos() + Float3{ 0, 0, 1 });
	cameraMain.setOrgPos(Trains[USER].carPos[atCar]);

	while (System::Update())
	{
		//メインカメラ描画
		{
			const ScopedRenderTarget3D rtmain{ rtexMain.clear(BGCOLOR) };

			Graphics3D::SetCameraTransform(cameraMain.getViewProj(), cameraMain.getEyePos());
			{
				const ScopedRenderStates3D rs{ SamplerState::RepeatAniso, RasterizerState::SolidCullBack };

				Trains[USER].progressT = AddLooped(Trains[USER].progressT, 1.0, Trains[USER].speedT);

				//曲線座標系
				for (int i = 0; i < Railways.size(); i+=3)
					updateTrain(Railways[i], Trains[i].progressT, Trains[i].spacing, Trains[i].carPos, Trains[i].carRotQ);

				//電車描画
				for (int i = 0; i < Railways.size(); i+=3)
					drawTrain(meshMaple, Trains[i].carPos, Trains[i].carRotQ);

				//線路描画
				drawMaple(meshMaple, mapSubway, chunkSubway.mat(),
						   Trains[USER].progressT, Railways[USER], Rect{ Arg::center(3, 0), 7, 80 });

				//カメラ制御
				if (atCar == CAR0) mapInterior.grid[13] = '#';	//先頭車両と後尾車両は通路を塞ぐ
				else if (atCar == CAR4) mapInterior.grid[400] = '#';
				else { mapInterior.grid[13] = ':'; mapInterior.grid[400] = ':'; }

				OrientedBox ob = meshMaple.getCollider(TRAIN_MID);
				updateMainCamera(cameraMain, Trains[USER].carPos[atCar], Trains[USER].carRotQ[atCar], mapInterior, ob);
			}

			RENDERTEXTURE(rtexMain, MAINWINDOW);

		}


#define TEST_SUBCAMERA
#ifdef TEST_SUBCAMERA

		static double altitude = 10;

		//サブカメラ描画
		{
#define TEST_BIRD
#ifdef TEST_BIRD
			cameraSub.setOrgPos(Trains[USER].carPos[atCar])
				.setEyePos(cameraMain.getEyePos().movedBy(-altitude, altitude, 0))
				.setFocusPos(cameraMain.getEyePos())
				.updateViewWorld(Trains[USER].carPos[atCar], Trains[USER].carRotQ[atCar])
				.updateViewProj();
#endif

			Graphics3D::SetCameraTransform(cameraSub.getViewProj(), cameraSub.getEyePos());
			{
				TRANSPARENTBG(rtexSub);

				//電車描画
				drawTrain(meshMaple, Trains[USER].carPos, Trains[USER].carRotQ);
				testLine(Railways[3], cameraMain, meshMapleSub);

				PixieMesh::drawGridBox();
			}

			RENDERTEXTURE(rtexSub, SUBWINDOW);
		}
#endif

#define USE_DEBUGGUI
#ifdef USE_DEBUGGUI
		{
			auto& c = cameraMain;
			double wheel = Mouse::Wheel() / 10;
			Rect base1 = Rect{ SWPOS.x,SWPOS.y,300,22 };
			{static double p; Trains[USER].speedT = guiValue(*fontS, U"speedT:{:.8f}"_fmt(Trains[USER].speedT), base1.moveBy(0, 20), wheel / 1000000, 0.0, Trains[USER].speedT, p, 0.0, 1.0); }
			{static double p; Trains[USER].progressT = guiValue(*fontS, U"progressT:{:.2f}"_fmt(Trains[USER].progressT), base1.moveBy(0, 20), wheel / 10, 0.0, Trains[USER].progressT, p, 0.0, 1.0); }
			{static int p; atCar = guiValue(*fontS, U"CAR[{}]"_fmt(atCar), base1.moveBy(0, 20), wheel, 0, atCar, p, 0, 4); }
			{static double p; altitude = guiValue(*fontS, U"Altitude:{:.1f}m"_fmt(altitude), base1.moveBy(0, 20), wheel*100, 0.0, altitude, p, 3.0, 10000.0); }

			{static Float3 p; c.m_orgLocal = guiValue(*fontS, U" org:{:.1f}"_fmt(c.m_orgLocal), base1.moveBy(0, 20), c.m_orgLocal, p); }
			{static Float3 p; c.m_eyeLocal = guiValue(*fontS, U" eyeL:{:.1f}"_fmt(c.m_eyeLocal), base1.moveBy(0, 20), c.m_eyeLocal, p); }
			{static Float3 p; c.m_focusLocal = guiValue(*fontS, U" focusL:{:.1f}"_fmt(c.m_focusLocal), base1.moveBy(0, 20), c.m_focusLocal, p); }
			{static Float3 p;                  guiValue(*fontS, U" eyeW:{:.1f}"_fmt(c.m_eyeWorld.xyz()), base1.moveBy(0, 20), c.m_eyeWorld.xyz(), p); }
			{static Float3 p;                  guiValue(*fontS, U" focusW:{:.1f}"_fmt(c.m_focusWorld.xyz()), base1.moveBy(0, 20), c.m_focusWorld.xyz(), p); }
			{static Float3 p;                  guiValue(*fontS, U"eyeT:{:.1f}"_fmt(c.m_orgLocal + c.m_eyeLocal), base1.moveBy(0, 20), c.m_orgLocal + c.m_eyeLocal, p); }
			{static Float3 p; Trains[USER].carPos[0] = guiValue(*fontS, U"CAR0:{:.1f}"_fmt(Trains[USER].carPos[0]), base1.moveBy(0, 20), Trains[USER].carPos[0], p); }
			{static Float3 p; Trains[USER].carPos[1] = guiValue(*fontS, U"CAR1:{:.1f}"_fmt(Trains[USER].carPos[1]), base1.moveBy(0, 20), Trains[USER].carPos[1], p); }
			{static Float3 p; Trains[USER].carPos[2] = guiValue(*fontS, U"CAR2:{:.1f}"_fmt(Trains[USER].carPos[2]), base1.moveBy(0, 20), Trains[USER].carPos[2], p); }
			{static Float3 p; Trains[USER].carPos[3] = guiValue(*fontS, U"CAR3:{:.1f}"_fmt(Trains[USER].carPos[3]), base1.moveBy(0, 20), Trains[USER].carPos[3], p); }
			{static Float3 p; Trains[USER].carPos[4] = guiValue(*fontS, U"CAR4:{:.1f}"_fmt(Trains[USER].carPos[4]), base1.moveBy(0, 20), Trains[USER].carPos[4], p); }
			{static Float3 p; Trains[USER].carPos[5] = guiValue(*fontS, U"CAR5:{:.1f}"_fmt(Trains[USER].carPos[5]), base1.moveBy(0, 20), Trains[USER].carPos[5], p); }

			Rect base2 = Rect{ SWPOS.x + 500,SWPOS.y,200,22 };
			OrientedBox ob = meshMaple.getCollider(COLLIDER_MID);
			Float2 newpos = cameraMain.getEyePos().xz() + Float2{ ob.w / 2, 0 };
			Float3 gridpos = {};
			mapInterior.GetChar(newpos, ob, gridpos);

			{static OrientedBox::position_type p; guiValue(*fontS, U"OB.size:{:.1f}"_fmt(ob.size), base2.moveBy(0, 20), ob.size, p); }
			{static Size p;   guiValue(*fontS, U"Interior.size:{}"_fmt(mapInterior.size), base2.moveBy(0, 20), mapInterior.size, p); }
			{static Float2 p; guiValue(*fontS, U"newpos:{:.1f}"_fmt(newpos), base2.moveBy(0, 20), newpos, p); }
			{static Float3 p; guiValue(*fontS, U"gridpos:{:.0f}"_fmt(gridpos), base2.moveBy(0, 20), gridpos, p); }

			SimpleGUI::VerticalSlider(Trains[USER].speedT, 0.0, 0.00010000, Vec2{ SWPOS.x + 750,SWPOS.y }, 200);

			Rect base3 = Rect{ SWPOS.x + 700,SWPOS.y-80 ,100,700 };

			int idx = gridpos.y * (mapInterior.size.x + 1) + gridpos.x;
			char buf[] = " ";
			buf[0] = mapInterior1.at(idx);
			mapInterior1.at(idx) = 'A';

			(*fontXS)(Unicode::Widen(mapInterior1)).draw(base3.moveBy(0, 20));

			mapInterior1.at(idx) = buf[0];

		}
#endif

	}
}

