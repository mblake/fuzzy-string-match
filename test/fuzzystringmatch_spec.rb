#!/usr/bin/env ruby
# -*- encoding: utf-8 -*-
#
# fuzzystringmatch_spec.rb -  "RSpec file for FuzzyStringMatch module "
#  
#   Copyright (c) 2010  Kiyoka Nishiyama  <kiyoka@sumibi.org>
#
#  Licensed to the Apache Software Foundation (ASF) under one or more
#  contributor license agreements.  See the NOTICE file distributed with
#  this work for additional information regarding copyright ownership.
#  The ASF licenses this file to You under the Apache License, Version 2.0
#  (the "License"); you may not use this file except in compliance with
#  the License.  You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
#
#  
require 'fuzzystringmatch'
require 'amatch'

describe FuzzyStringMatch, "when some string distances (Pure) are" do
  before do
    @jarow = FuzzyStringMatch::JaroWinkler.new.create
  end
  it "should" do
    @jarow.getDistance( "henka",     "henkan"    ).should be_within(0.0001).of(0.9722)
    @jarow.getDistance( "al",        "al"        ).should == 1.0
    @jarow.getDistance( "martha",    "marhta"    ).should be_within(0.0001).of(0.9611)
    @jarow.getDistance( "jones",     "johnson"   ).should be_within(0.0001).of(0.8323)
    @jarow.getDistance( "abcvwxyz",  "cabvwxyz"  ).should be_within(0.0001).of(0.9583)
    @jarow.getDistance( "dwayne",    "duane"     ).should be_within(0.0001).of(0.8400)
    @jarow.getDistance( "dixon",     "dicksonx"  ).should be_within(0.0001).of(0.8133)
    @jarow.getDistance( "fvie",      "ten"       ).should == 0.0
    lambda {
      d1 = @jarow.getDistance("zac ephron", "zac efron")
      d2 = @jarow.getDistance("zac ephron", "kai ephron")
      d1 > d2
    }.should be_true
    lambda {
      d1 = @jarow.getDistance("brittney spears", "britney spears")
      d2 = @jarow.getDistance("brittney spears", "brittney startzman")
      d1 > d2
    }.should be_true
  end
end

describe FuzzyStringMatch, "when some string distances (Native) are" do
  before do
    @jarow = FuzzyStringMatch::JaroWinkler.new.create( :native )
  end
  it "should" do
    @jarow.getDistance( "henka",     "henkan"    ).should be_within(0.0001).of(0.9722)
    @jarow.getDistance( "al",        "al"        ).should == 1.0
    @jarow.getDistance( "martha",    "marhta"    ).should be_within(0.0001).of(0.9611)
    @jarow.getDistance( "jones",     "johnson"   ).should be_within(0.0001).of(0.8323)
    @jarow.getDistance( "abcvwxyz",  "cabvwxyz"  ).should be_within(0.0001).of(0.9583)
    @jarow.getDistance( "dwayne",    "duane"     ).should be_within(0.0001).of(0.8400)
    @jarow.getDistance( "dixon",     "dicksonx"  ).should be_within(0.0001).of(0.8133)
    @jarow.getDistance( "fvie",      "ten"       ).should == 0.0
    lambda {
      d1 = @jarow.getDistance("zac ephron", "zac efron")
      d2 = @jarow.getDistance("zac ephron", "kai ephron")
      d1 > d2
    }.should be_true
    lambda {
      d1 = @jarow.getDistance("brittney spears", "britney spears")
      d2 = @jarow.getDistance("brittney spears", "brittney startzman")
      d1 > d2
    }.should be_true
  end
end

def amatch_getDistance( s1, s2 )
  @jarow = Amatch::JaroWinkler.new( s1 )
  @jarow.match( s2 )
end


describe FuzzyStringMatch, "when some UTF8 string distances (Pure) are" do
  before do
    @jarow = FuzzyStringMatch::JaroWinkler.new.create
  end
  it "should" do
    @jarow.getDistance( "ａｌ",              "ａｌ"             ).should == 1.0
    @jarow.getDistance( "ｍａｒｔｈａ",      "ｍａｒｈｔａ"     ).should be_within(0.0001).of(0.9611)
    @jarow.getDistance( "ｊｏｎｅｓ",        "ｊｏｈｎｓｏｎ"   ).should be_within(0.0001).of(0.8323)
    @jarow.getDistance( "ａｂｃｖｗｘｙｚ",  "ｃａｂｖｗｘｙｚ" ).should be_within(0.0001).of(0.9583)
    @jarow.getDistance( "ｄｗａｙｎｅ",      "ｄｕａｎｅ"       ).should be_within(0.0001).of(0.8400)
    @jarow.getDistance( "ｄｉｘｏｎ",        "ｄｉｃｋｓｏｎｘ" ).should be_within(0.0001).of(0.8133)
    @jarow.getDistance( "ｆｖｉｅ",          "ｔｅｎ"           ).should == 0.0
    lambda {
      d1 = @jarow.getDistance("ｚａｃ　ｅｐｈｒｏｎ", "ｚａｃ　ｅｆｒｏｎ")
      d2 = @jarow.getDistance("ｚａｃ　ｅｐｈｒｏｎ", "ｋａｉ　ｅｐｈｒｏｎ")
      d1 > d2
    }.should be_true
    lambda {
      d1 = @jarow.getDistance("ｂｒｉｔｔｎｅｙ　ｓｐｅａｒｓ", "ｂｒｉｔｎｅｙ　ｓｐｅａｒｓ")
      d2 = @jarow.getDistance("ｂｒｉｔｔｎｅｙ　ｓｐｅａｒｓ", "ｂｒｉｔｔｎｅｙ　ｓｔａｒｔｚｍａｎ")
      d1 > d2
    }.should be_true
    @jarow.getDistance( "スパゲティー",              "スパゲッティー"             ).should be_within(0.0001).of(0.9666)
    @jarow.getDistance( "スパゲティー",              "スパゲティ"                 ).should be_within(0.0001).of(0.9722)
    @jarow.getDistance( "スティービー・ワンダー",    "スピーディー・ワンダー"     ).should be_within(0.0001).of(0.8561)
    @jarow.getDistance( "マイケル・ジャクソン",      "ジャイケル・マクソン"       ).should be_within(0.0001).of(0.8000)
    @jarow.getDistance( "まつもとゆきひろ",          "まつもとひろゆき"           ).should be_within(0.0001).of(0.9500)
    @jarow.getDistance( "クライエント",              "クライアント"               ).should be_within(0.0001).of(0.9222)
    @jarow.getDistance( "サーバー",                  "サーバ"                     ).should be_within(0.0001).of(0.9416)
  end
end


describe Amatch, "when use Amatch gem, results are" do
  it "should" do
    amatch_getDistance( "henka",     "henkan"    ).should be_within(0.0001).of(0.9666)  ## amatch's result value is different from lucene version.
    amatch_getDistance( "al",        "al"        ).should == 1.0
    amatch_getDistance( "martha",    "marhta"    ).should be_within(0.0001).of(0.9611)
    amatch_getDistance( "jones",     "johnson"   ).should be_within(0.0001).of(0.8323)
    amatch_getDistance( "abcvwxyz",  "cabvwxyz"  ).should be_within(0.0001).of(0.9583)
    amatch_getDistance( "dwayne",    "duane"     ).should be_within(0.0001).of(0.8400)
    amatch_getDistance( "dixon",     "dicksonx"  ).should be_within(0.0001).of(0.8133)
    amatch_getDistance( "fvie",      "ten"       ).should == 0.0
    lambda {
      d1 = amatch_getDistance("zac ephron", "zac efron")
      d2 = amatch_getDistance("zac ephron", "kai ephron")
      d1 > d2
    }.should be_true
    lambda {
      d1 = amatch_getDistance("brittney spears", "britney spears")
      d2 = amatch_getDistance("brittney spears", "brittney startzman")
      d1 > d2
    }.should be_true
  end
end
