# include <Siv3D.hpp>   

typedef uint32 uint;

# include "PixieMesh.hpp"
# include "PixieCamera.hpp"
# include "LineString3D.hpp"

# include "Maple2.hpp"
using namespace libMaple;

constexpr ColorF        GRAY1 = { 0.01, 0.01, 0.01, 0.5 };
constexpr ColorF        GRAY3 = { 0.03, 0.03, 0.03, 0 };
constexpr ColorF        GRAY5 = { 0.05, 0.05, 0.05, 0 };
constexpr ColorF        GRAY8 = { 0.08, 0.08, 0.08, 0 };

constexpr ColorF        BGCOLOR = GRAY5;
constexpr TextureFormat TEXFORMAT = TextureFormat::R8G8B8A8_Unorm_SRGB;
constexpr Size			WINDOWSIZE = { 1920, 768 };

constexpr StringView	APATH = U"../Asset/";

#define TRANSPARENTBG(RTEX)	BlendState blend = BlendState::Default2D;\
							blend.srcAlpha = Blend::SrcAlpha;\
							blend.dstAlpha = Blend::DestAlpha;\
							blend.opAlpha = BlendOp::Max;\
							const ScopedRenderStates3D bs{ BlendState::OpaqueAlphaToCoverage,\
								  SamplerState::RepeatAniso, RasterizerState::SolidCullNone };\
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
	return value;
}
template <typename V>
inline V guiValue(const Font& font, const String& text, Rect rect, V value, V& prev)
{
	font(text).draw(rect, rect.mouseOver() ? ColorF(Palette::Red) : ColorF(Palette::White));
	return value;
}

template <typename V>
inline V AddLooped(V value, const V& max,const V& amount)
{
	value += amount;
	if (value > max) value -= max;
	if (value < 0.0) value += max;
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
	uint32 shiftTime = 10000;
	double changeTime = 0.000000;
	double dstSpeedT = 0.000001;
	double srcSpeedT = 0.000000;

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


enum {
	TRAIN_TAIL = 29, TRAIN_MID = 28, TRAIN_HEAD = 27,
	COLLIDER_TAIL = 32, COLLIDER_MID = 31, COLLIDER_HEAD = 30,
	ARROW = 33, CAMERA = 34,
};

Array<Entity> entMaple = {
	Entity{ "A1",  0}, Entity{ "D6", 1 }, Entity{ "A2", 2 }, Entity{ "A3", 3 },
	Entity{ "B1",  4}, Entity{ "B2", 5 }, Entity{ "B3", 6 }, Entity{ "C1", 7 },
	Entity{ "C2",  8}, Entity{ "C3", 9 }, Entity{ "C4", 10}, Entity{ "C5", 11},
	Entity{ "C6", 12}, Entity{ "D1", 13}, Entity{ "D2", 14}, Entity{ "D3", 19},
	Entity{ "D4", 15}, Entity{ "D5", 16}, Entity{ "ES", 17},

	Entity{ "SW", 18 },
	Entity{ "R1", 20 }, Entity{ "R2", 21}, Entity{ "R3", 22},
	Entity{ "S1", 23 }, Entity{ "S2", 24}, Entity{ "S3", 25},

	Entity{ "DR", 26}, Entity{ "TR1", TRAIN_HEAD },Entity{ "TR2", TRAIN_MID },Entity{ "TR3", TRAIN_TAIL },

	Entity{ "CL1", COLLIDER_TAIL}, Entity{ "CL2", COLLIDER_MID}, Entity{ "CL3", COLLIDER_HEAD},

	Entity{ "AR", ARROW},	Entity{ "CAM", CAMERA},
};

void updateMainCamera(PixieCamera& camera, const Float3& orgpos, const Quaternion& orgrotq,
					   const Maple& maple, const OrientedBox& ob)
{
	Float2 delta = Cursor::DeltaF();
	Float2 point2D = Cursor::PosF();
	Float3 distance = {};

	float speedM = 1;
	if (MouseM.down())								
	{
		const Ray mouseRay = camera.screenToRay(Cursor::PosF());
		if (mouseRay.direction.hasNaN()) assert(1);
	}

	if (MouseM.pressed())
	{
		camera.panX(delta.x * speedM / 100);
		camera.panY(delta.y * -speedM / 100);
	}

	if (MouseR.pressed())							
	{
		camera.trackX(delta.x * speedM / 20);
		camera.craneY(delta.y * speedM / 20);
	}

	if (KeyLControl.pressed() || MouseM.pressed())	
	{
		Float3 eyepos = cameraMain.getOrgEyePos();
		Float3 focuspos = cameraMain.getFocusPos();
		bool rev = camera.dolly((float)Mouse::Wheel() * speedM / 5, true);	
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

	mesh.setMove(camera.getEyeWorld().xyz() )
		.setRotateQ(camera.getQForward() * Trains[USER].carRotQ[atCar])
		.drawMesh(ARROW, Palette::Green);
}


void updateTrain( Array<LineString3D>& railways, Array<Train> &trains )
{
	for (int t = 0; t < trains.size(); t++)
	{
		Train& TR = trains[t];
		TR.progressT = AddLooped(TR.progressT, 1.0, TR.speedT);

		for (int i = 0; i < TR.carPos.size(); i++)
			TR.carPos[i] = railways[t].getCurvedPoint(AddLooped(TR.progressT, 1.0, TR.spacing * (i + 2)));

		for (int i = 0; i < TR.carRotQ.size() - 1; i++)
			TR.carRotQ[i] = PixieCamera::getQLookAt(TR.carPos[i], TR.carPos[i + 1]);

		if (t != USER)
		{
			if (TR.shiftTime)
			{
				if (--TR.shiftTime == 0)
				{
					TR.changeTime = 1.0;
					TR.dstSpeedT = (Random(1.0) - 0.5) / 10000;	
					TR.srcSpeedT = TR.speedT;
				}
			}

			if (TR.changeTime > 0.0)
			{
				TR.speedT = Math::Lerp(TR.dstSpeedT, TR.srcSpeedT, TR.changeTime);

				TR.changeTime -= 0.01;
				if (TR.changeTime <= 0.0)
				{
					TR.shiftTime = Random(500) + 500;
					TR.speedT = TR.dstSpeedT;		
				}
			}
		}
	}

}

void changeCar( int car, int meshid, PixieMesh& meshMaple, const Array<Float3>& carPos,
	          const Array<Quaternion>& carRotQ)
{
	OrientedBox ob = meshMaple.getCollider(meshid);
	ob = Geometry3D::TransformBoundingOrientedBox(ob, meshMaple.getMat());
	if (ob.contains(cameraMain.getEyeWorld().xyz()))
	{
		if (atCar != car)
		{
			Float3 dir = cameraMain.getFocusPos() - cameraMain.getEyePos();

			Float3 eyepos = Mat4x4(carRotQ[car]).inverse()
				.transformPoint(cameraMain.getEyeWorld().xyz() - carPos[car]);

			cameraMain.setEyePos(eyepos);
			cameraMain.setFocusPos(eyepos + dir);
			atCar = car;
		}
	}

	if (atCar == car) meshMapleSub = meshMaple;
}


void drawTrain(PixieMesh& meshMaple, int t, Array<Train>& trains )
{
	const Array<int> car{ TRAIN_TAIL, TRAIN_MID, TRAIN_MID, TRAIN_MID, TRAIN_HEAD };

	Train& TR = trains[t];

	const double& pr = trains[USER].progressT ;

	for (int i = 0; i < car.size(); i++)
	{
			
		meshMaple.setMove(TR.carPos[i]).setRotateQ(TR.carRotQ[i]).drawMesh(car[i]);

		if (t == USER)
		{
			changeCar(i, car[i], meshMaple, TR.carPos, TR.carRotQ);
			cameraMain.setOrgPos(TR.carPos[atCar]);			 
		}
	}
}

void Main()
{
	const Font FONTS(15); fontS = &FONTS;
	const Font FONTXS(8, FileSystem::GetFolderPath(SpecialFolder::SystemFonts) + U"MSGOTHIC.TTC"); fontXS = &FONTXS;
	static MSRenderTexture rtexMain = { (unsigned)WINDOWSIZE.x, (unsigned)WINDOWSIZE.y, TEXFORMAT, HasDepth::Yes };
	static MSRenderTexture rtexSub = { (unsigned)WINDOWSIZE.x, (unsigned)WINDOWSIZE.y, TEXFORMAT, HasDepth::Yes };

	Window::Resize(WINDOWSIZE);
	Window::SetStyle(WindowStyle::Sizable);

	const Size SWPOS = Window::GetState().virtualSize - Point{ 805,365 };
	const RectF SUBWINDOW = RectF{ SWPOS, 800,360 };
	const RectF MAINWINDOW = RectF{ 0,0, WINDOWSIZE };

	constexpr double RADIUS = 1000;
	constexpr double VOLATILITY = 500;

	const Array<std::string> mapSubway1 {
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

	std::string mapInterior1 =
	{        
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

	LineString3D curvedway = trackway.catmullRomClosed(80);						
	double len = curvedway.updateLengthList();

	LineString3D railway{ curvedway.getEquatePoints((int)(len / 4 )) };			
	railway.updateLengthList();

	Array<LineString3D> Railways(NUM_TRAINS);									
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

	for (int i = 0; i < Railways.size(); i++)		
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

		Trains[i].shiftTime = Random(500) + 500;	
		Trains[i].speedT = (Random(1.0)-0.5) / 1000;
		Trains[i].progressT = Random(1.0) ;			
	}

	std::string subwaymap = {};

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

	PixieMesh meshMaple; meshMaple.initModel(APATH + U"Subway.07.glb", MODELNOA, WINDOWSIZE, USE_INDIVIDUAL);
	PixieMesh meshFont;	meshFont.initModel(APATH + U"Font.01.glb", MODELNOA, WINDOWSIZE, USE_INDIVIDUAL);

	meshMaple.obCollider[TRAIN_HEAD].size += { 0, 0, -1.0};
	meshMaple.obCollider[TRAIN_MID].size += { 0, 0, -1.0};
	meshMaple.obCollider[TRAIN_TAIL].size += { 0, 0, -1.0};

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
		{
			const ScopedRenderTarget3D rtmain{ rtexMain.clear(BGCOLOR) };

			Graphics3D::SetCameraTransform(cameraMain.getViewProj(), cameraMain.getEyePos());
			{
				const ScopedRenderStates3D rs{ SamplerState::RepeatAniso, RasterizerState::SolidCullNone };

				updateTrain(Railways, Trains);

				for (int i = 0; i < Railways.size(); i += 1)
					drawTrain(meshMaple, i, Trains );

				drawMaple(meshMaple, mapSubway, chunkSubway.mat(),
					Trains[USER].progressT, Railways[USER], Rect{ Arg::center(3, 0), 7, 80 });

				static float AAA = 0.0;
				meshFont.setCamera(cameraMain);
				meshFont.setMove( cameraMain.getEyeWorld().xyz() ).setScale({ 2,2,2 })
					.setRotateQ(PixieCamera::getQLookAt( cameraMain.getFocusWorld() , cameraMain.getEyeWorld().xyz()) )
					.drawString(U"ABCDEFGHIJKLMNOPQRSTUVWXYZ\n"
						         "abcdefghijklmnopqrstuvwxyz\n"
						"0123456789!\"#$%&'()=~|[]`*{}", 1.5, 4.0, { 0.9,0.9,0.9 }, 0, 83);
				AAA += 0.1;
				if (AAA > 83.0) AAA = 0.0;

				if (atCar == CAR0) mapInterior.grid[13] = '#';	
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

				meshFont.drawString(U"ABCDEFGHIJKLMNOPQRSTUVWXYZ\n"
						"abcdefghijklmnopqrstuvwxyz\n"
						"0123456789!\"#$%&'()=~|[]`*{}", 1.5, 4.0, { 0.9,0.9,0.9 }, 0, 83);

				drawTrain(meshMaple, USER, Trains);
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

			Rect base2 = Rect{ SWPOS.x + 250,SWPOS.y,200,22 };
			OrientedBox ob = meshMaple.getCollider(COLLIDER_MID);
			Float2 newpos = cameraMain.getEyePos().xz() + Float2{ ob.w / 2, 0 };
			Float3 gridpos = {};
			mapInterior.GetChar(newpos, ob, gridpos);

			{static OrientedBox::position_type p; guiValue(*fontS, U"OB.size:{:.1f}"_fmt(ob.size), base2.moveBy(0, 20), ob.size, p); }
			{static Size p;   guiValue(*fontS, U"Interior.size:{}"_fmt(mapInterior.size), base2.moveBy(0, 20), mapInterior.size, p); }
			{static Float2 p; guiValue(*fontS, U"newpos:{:.1f}"_fmt(newpos), base2.moveBy(0, 20), newpos, p); }
			{static Float3 p; guiValue(*fontS, U"gridpos:{:.0f}"_fmt(gridpos), base2.moveBy(0, 20), gridpos, p); }

			{static double p; Trains[0].progressT = guiValue(*fontS, U"Train0:{:.4f}"_fmt(Trains[0].progressT), base2.moveBy(0, 20), wheel / 10, 0.0, Trains[0].progressT, p, 0.0, 1.0); }
			{static double p; Trains[1].progressT = guiValue(*fontS, U"Train1:{:.4f}"_fmt(Trains[1].progressT), base2.moveBy(0, 20), wheel / 10, 0.0, Trains[1].progressT, p, 0.0, 1.0); }
			{static double p; Trains[2].progressT = guiValue(*fontS, U"Train2:{:.4f}"_fmt(Trains[2].progressT), base2.moveBy(0, 20), wheel / 10, 0.0, Trains[2].progressT, p, 0.0, 1.0); }
			{static double p; Trains[4].progressT = guiValue(*fontS, U"Train4:{:.4f}"_fmt(Trains[4].progressT), base2.moveBy(0, 20), wheel / 10, 0.0, Trains[4].progressT, p, 0.0, 1.0); }
			{static double p; Trains[5].progressT = guiValue(*fontS, U"Train5:{:.4f}"_fmt(Trains[5].progressT), base2.moveBy(0, 20), wheel / 10, 0.0, Trains[5].progressT, p, 0.0, 1.0); }
			{static double p; Trains[6].progressT = guiValue(*fontS, U"Train6:{:.4f}"_fmt(Trains[6].progressT), base2.moveBy(0, 20), wheel / 10, 0.0, Trains[6].progressT, p, 0.0, 1.0); }

			{static double p; Trains[0].speedT = guiValue(*fontS, U"Speed0:{:.9f}"_fmt(Trains[0].speedT), base2.moveBy(0, 20), wheel / 1000000, 0.0, Trains[0].speedT, p, 0.0, 1.0); }
			{static double p; Trains[1].speedT = guiValue(*fontS, U"Speed1:{:.9f}"_fmt(Trains[1].speedT), base2.moveBy(0, 20), wheel / 1000000, 0.0, Trains[1].speedT, p, 0.0, 1.0); }
			{static double p; Trains[2].speedT = guiValue(*fontS, U"Speed2:{:.9f}"_fmt(Trains[2].speedT), base2.moveBy(0, 20), wheel / 1000000, 0.0, Trains[2].speedT, p, 0.0, 1.0); }
			{static double p; Trains[4].speedT = guiValue(*fontS, U"Speed4:{:.9f}"_fmt(Trains[4].speedT), base2.moveBy(0, 20), wheel / 1000000, 0.0, Trains[4].speedT, p, 0.0, 1.0); }
			{static double p; Trains[5].speedT = guiValue(*fontS, U"Speed5:{:.9f}"_fmt(Trains[5].speedT), base2.moveBy(0, 20), wheel / 1000000, 0.0, Trains[5].speedT, p, 0.0, 1.0); }
			{static double p; Trains[6].speedT = guiValue(*fontS, U"Speed6:{:.9f}"_fmt(Trains[6].speedT), base2.moveBy(0, 20), wheel / 1000000, 0.0, Trains[6].speedT, p, 0.0, 1.0); }

			Rect base3 = Rect{ SWPOS.x + 450,SWPOS.y,200,22 };
			{static uint32 p; Trains[0].shiftTime = guiValue(*fontS, U"shiftT0:{}"_fmt(Trains[0].shiftTime), base3.moveBy(0, 20), Trains[0].shiftTime, p); }
			{static uint32 p; Trains[1].shiftTime = guiValue(*fontS, U"shiftT1:{}"_fmt(Trains[1].shiftTime), base3.moveBy(0, 20), Trains[1].shiftTime, p); }
			{static uint32 p; Trains[2].shiftTime = guiValue(*fontS, U"shiftT2:{}"_fmt(Trains[2].shiftTime), base3.moveBy(0, 20), Trains[2].shiftTime, p); }
			{static uint32 p; Trains[4].shiftTime = guiValue(*fontS, U"shiftT4:{}"_fmt(Trains[4].shiftTime), base3.moveBy(0, 20), Trains[4].shiftTime, p); }
			{static uint32 p; Trains[5].shiftTime = guiValue(*fontS, U"shiftT5:{}"_fmt(Trains[5].shiftTime), base3.moveBy(0, 20), Trains[5].shiftTime, p); }
			{static uint32 p; Trains[6].shiftTime = guiValue(*fontS, U"shiftT6:{}"_fmt(Trains[6].shiftTime), base3.moveBy(0, 20), Trains[6].shiftTime, p); }

			{static double p; Trains[0].changeTime = guiValue(*fontS, U"changeT0:{:.3f}"_fmt(Trains[0].changeTime), base3.moveBy(0, 20), Trains[0].changeTime, p); }
			{static double p; Trains[1].changeTime = guiValue(*fontS, U"changeT1:{:.3f}"_fmt(Trains[1].changeTime), base3.moveBy(0, 20), Trains[1].changeTime, p); }
			{static double p; Trains[2].changeTime = guiValue(*fontS, U"changeT2:{:.3f}"_fmt(Trains[2].changeTime), base3.moveBy(0, 20), Trains[2].changeTime, p); }
			{static double p; Trains[4].changeTime = guiValue(*fontS, U"changeT4:{:.3f}"_fmt(Trains[4].changeTime), base3.moveBy(0, 20), Trains[4].changeTime, p); }
			{static double p; Trains[5].changeTime = guiValue(*fontS, U"changeT5:{:.3f}"_fmt(Trains[5].changeTime), base3.moveBy(0, 20), Trains[5].changeTime, p); }
			{static double p; Trains[6].changeTime = guiValue(*fontS, U"changeT6:{:.3f}"_fmt(Trains[6].changeTime), base3.moveBy(0, 20), Trains[6].changeTime, p); }

			Rect base4 = Rect{ SWPOS.x + 700,SWPOS.y - 80 ,100,700 };
			(*fontXS)(Unicode::Widen(mapInterior1)).draw(base4.moveBy(0, 20));

			SimpleGUI::VerticalSlider(Trains[USER].speedT, 0.0, 0.00010000, Vec2{ SWPOS.x + 750,SWPOS.y }, 200);
			int idx = gridpos.y * (mapInterior.size.x + 1) + gridpos.x;
			char buf[] = " ";
			buf[0] = mapInterior1.at(idx);
			mapInterior1.at(idx) = 'A';



			mapInterior1.at(idx) = buf[0];

		}
#endif

	}
}

